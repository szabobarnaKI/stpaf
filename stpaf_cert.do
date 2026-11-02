* stpaf_cert.do  --  certification script for stpaf 0.4.0
* Run with stpaf.ado and stpaf.sthlp on the adopath (e.g. -adopath ++ "path/to/folder"-).
* NOTE: written without access to a Stata session; please run it once and report any failure.
*
* What is checked
*   1. E (model-expected events) equals the observed number of failures when nothing is restricted
*   2. A scenario that does not change the linear predictor gives PAF = 0 and SE = 0
*   3. Analytic case: eliminating the only exposure within exposed person-time gives 1 - exp(-b)
*   4. set(z=0) and eliminate(z) agree; each reproduces the joint scenario when only one change is given
*   5. IF-based SE agrees with a clustered bootstrap (recurrent events, time-varying exposure, strata)
*   6. Horizon results are consistent: cumulative PAF at the maximum time equals the default
*   7. fweights reproduce the expanded data set; unit pweights reproduce the unweighted result
*   8. IF-based SE agrees with the cluster jackknife (deterministic check of the variance code)
*   9. vce(jackknife) refit reproduces the explicit leave-one-cluster-out loop of test 8;
*      the Newton jackknife (default) equals the refits up to stcox's convergence tolerance,
*      also for several scenarios and horizons; warm and cold starts (nowarm) agree to 1e-9
*  10. the same without strata (0.3.1 fix: vce(jackknife) failed for unstratified models)
*  11. cluster(): precedence over e(clustvar)/stset id, coarser and non-nested clusters,
*      missing values, and data stset like ghoshlin output ([pw=weight_c] without id())
*  12. ties: efron without tied failure times uses the Newton jackknife; efron with tied
*      failure times falls back to stcox refits; breslow with ties: Newton = refits

version 15.0
clear all
set more off
set type double
set seed 20261009

capture program drop _cert_sim
program define _cert_sim
    syntax, n(integer)
    clear
    set obs `n'
    gen long id = _n
    gen double age = rnormal()
    gen byte stratum = 1 + (runiform() > .5)
    gen double C = 2 + 3*runiform()
    gen double sw = rexponential(2)
    gen double fr = rgamma(2, .5)
    gen double len0 = min(sw, C)
    gen double len1 = max(C - sw, 0)
    gen double r0 = (0.6 + 0.4*(stratum == 2))*fr*exp(0.3*age)
    gen double r1 = r0*exp(0.7)
    gen long n0 = rpoisson(r0*len0)
    gen long n1 = rpoisson(r1*len1)
    expand n0 + n1 + 1
    bysort id: gen long k = _n
    gen double tev = .
    replace tev = runiform()*len0 if k <= n0
    replace tev = sw + runiform()*len1 if k > n0 & k <= n0 + n1
    replace tev = C if k == n0 + n1 + 1
    gen byte ev = k <= n0 + n1
    sort id tev
    by id: gen double tstart = cond(_n == 1, 0, tev[_n-1])
    drop if tev <= tstart
    stset tev, id(id) failure(ev) time0(tstart) exit(time .)
    stsplit blk, after(time = sw) at(0)
    gen byte z = (_t0 >= sw)
    gen double t0b = _t0
    gen double t1b = _t
    gen byte db = _d
    // the simulation must give recurrent events (several failures per subject)
    tempvar nf
    bysort id: egen `nf' = total(_d)
    summarize `nf', meanonly
    assert r(max) > 1
end

capture program drop _cert_bs
program define _cert_bs, rclass
    // [if] is required by jackknife, which drops one cluster at a time via an if condition
    syntax [if]
    preserve
    if `"`if'"' != "" qui keep `if'
    stset t1b, id(newid) failure(db) time0(t0b) exit(time .)
    qui stcox z age, strata(stratum) vce(cluster newid) breslow
    qui stpaf, eliminate(z)
    return scalar paf = r(paf)
    restore
end

capture program drop _cert_bs0
program define _cert_bs0, rclass
    // as _cert_bs, without strata
    syntax [if]
    preserve
    if `"`if'"' != "" qui keep `if'
    stset t1b, id(newid) failure(db) time0(t0b) exit(time .)
    qui stcox z age, vce(cluster newid) breslow
    qui stpaf, eliminate(z)
    return scalar paf = r(paf)
    restore
end

_cert_sim, n(400)
qui count if _d
assert r(N) > 400
di as txt "recurrent-event data: " r(N) " failures in 400 subjects"
qui stcox z age, strata(stratum) vce(cluster id) breslow
local bz = _b[z]

* ---- 1. expected events equal observed failures
stpaf, eliminate(z)
assert reldif(r(E), r(N_fail)) < 1e-8
assert r(lb) < r(paf) & r(paf) < r(ub)
local paf1 = r(paf)
local se1  = r(se)
di as txt "[ok] 1: E = N_fail = " r(N_fail)

* ---- 2. null scenario
qui stcox age, strata(stratum) vce(cluster id) breslow
stpaf, eliminate(z)
assert abs(r(paf)) < 1e-10
assert r(se) < 1e-8
di as txt "[ok] 2: scenario outside the model gives PAF = 0"

* ---- 3. analytic case
qui stcox z age, strata(stratum) vce(cluster id) breslow
local bz = _b[z]
stpaf if z == 1, eliminate(z)
assert reldif(r(paf), 1 - exp(-`bz')) < 1e-8
di as txt "[ok] 3: PAF among exposed person-time = 1 - exp(-b) = " 1 - exp(-`bz')

* ---- 4. equivalent specifications
stpaf, eliminate(z)
local a = r(paf)
local sa = r(se)
stpaf, set(z=0)
assert reldif(r(paf), `a') < 1e-10 & reldif(r(se), `sa') < 1e-10
stpaf, eliminate(z) each
matrix T = r(table)
assert rowsof(T) == 1
stpaf, eliminate(z) cap(age=0) each
matrix T = r(table)
assert rowsof(T) == 3
assert reldif(T[2, 3], `a') < 1e-10
di as txt "[ok] 4: equivalent specifications agree"

* ---- 6. horizons
stpaf, eliminate(z) horizon(1 2 3)
matrix T = r(table)
assert rowsof(T) == 3
summ _t, meanonly
stpaf, eliminate(z) horizon(`=r(max)+1')
assert reldif(r(paf), `a') < 1e-10
di as txt "[ok] 6: horizon beyond the last event reproduces the default"

* ---- 5. bootstrap comparison
qui stcox z age, strata(stratum) vce(cluster id) breslow
timer clear
timer on 1
qui stpaf, eliminate(z)
timer off 1
local pafif = r(paf)
local seif  = r(se)
timer on 2
qui bootstrap paf = r(paf), reps(200) cluster(id) idcluster(newid) seed(7) nodots: _cert_bs
timer off 2
matrix V = e(V)
local sebs = sqrt(V[1, 1])
qui timer list
di as txt "IF  : PAF = " %7.4f `pafif' "  SE = " %7.4f `seif'  "  time = " %6.2f r(t1) "s"
di as txt "Boot: SE = " %7.4f `sebs' "  time = " %6.2f r(t2) "s (200 replications)"
assert inrange(`seif'/`sebs', 0.80, 1.25)
di as txt "[ok] 5: influence-function SE agrees with the clustered bootstrap"

* ---- 8. cluster jackknife (deterministic; close to the infinitesimal jackknife
*         that the influence function approximates, so it isolates coding errors
*         from finite-sample bootstrap behaviour)
*         written as an explicit leave-one-cluster-out loop (jackknife's own
*         e(sample) detection is confused by the preserve/restore in _cert_bs)
qui stcox z age, strata(stratum) vce(cluster id) breslow
qui stpaf, eliminate(z)
local seif = r(se)
qui summarize id, meanonly
local G = r(max)
capture drop newid
qui gen long newid = id
matrix JK = J(`G', 1, .)
timer on 3
forvalues c = 1/`G' {
    qui _cert_bs if id != `c'
    matrix JK[`c', 1] = r(paf)
}
timer off 3
drop newid
mata: st_local("sejk", strofreal(sqrt((rows(st_matrix("JK")) - 1)/rows(st_matrix("JK")) * ///
    sum((st_matrix("JK") :- mean(st_matrix("JK"))):^2)), "%21.15g"))
qui timer list
di as txt "Jackknife: SE = " %7.4f `sejk' "  time = " %6.2f r(t3) "s  (ratio IF/jackknife = " ///
    %5.3f `seif'/`sejk' ")"
assert inrange(`seif'/`sejk', 0.90, 1.05)
di as txt "[ok] 8: influence-function SE agrees with the cluster jackknife"

* ---- 9. vce(jackknife): refit = explicit loop of test 8; Newton = refit; warm = cold
qui stcox z age, strata(stratum) vce(cluster id) breslow
timer on 4
stpaf, eliminate(z) vce(jackknife) refit nodots
timer off 4
local serf = r(se)
assert "`r(jk_method)'" == "refit"
assert reldif(`serf', `sejk') < 1e-6
timer on 5
stpaf, eliminate(z) vce(jackknife)
timer off 5
local senw = r(se)
local slnw = r(se_logR)
local itw  = r(jk_iter)
assert "`r(jk_method)'" == "newton"
assert r(lb) < r(paf) & r(paf) < r(ub)
assert "`r(vce)'" == "jackknife"
// stcox refits stop at stcox's default convergence criteria (about 1e-6 relative)
assert reldif(`senw', `serf') < 1e-5
timer on 6
qui stpaf, eliminate(z) vce(jackknife) nowarm
timer off 6
assert "`r(jk_method)'" == "newton, nowarm"
assert reldif(r(se), `senw') < 1e-9 & reldif(r(se_logR), `slnw') < 1e-9
assert r(jk_iter) > `itw'
local itc = r(jk_iter)
* several scenarios and horizons, all columns of r(table)
qui stpaf, eliminate(z) cap(age=0) each horizon(1 2.5) vce(jackknife)
matrix TN = r(table)
qui stpaf, eliminate(z) cap(age=0) each horizon(1 2.5) vce(jackknife) nowarm
matrix TC = r(table)
qui stpaf, eliminate(z) cap(age=0) each horizon(1 2.5) vce(jackknife) refit nodots
matrix TR = r(table)
assert rowsof(TN) == 6
assert mreldif(TN, TC) < 1e-9
matrix SN = TN[1..., 4..6], TN[1..., 9]
matrix SR = TR[1..., 4..6], TR[1..., 9]
matrix PN = TN[1..., 3]
matrix PR = TR[1..., 3]
assert mreldif(PN, PR) < 1e-12            // point estimates use e(b) in both
assert mreldif(SN, SR) < 1e-5             // se, lb, ub, se_logR
qui timer list
di as txt "jackknife SE: explicit loop " %9.6f `sejk' "   refit " %9.6f `serf' "   Newton " %9.6f `senw'
di as txt "time: refit " %6.2f r(t4) "s   Newton (warm, " %4.2f `itw' " it.) " %6.2f r(t5) ///
    "s   Newton (nowarm, " %4.2f `itc' " it.) " %6.2f r(t6) "s"
di as txt "[ok] 9: refit = explicit loop; Newton = refit (to stcox's tolerance); warm = cold start"

* ---- 10. unstratified model: vce(jackknife) = explicit loop
qui stcox z age, vce(cluster id) breslow
qui stpaf, eliminate(z)
local paf0 = r(paf)
local seif0 = r(se)
capture drop newid
qui gen long newid = id
matrix JK0 = J(`G', 1, .)
forvalues c = 1/`G' {
    qui _cert_bs0 if id != `c'
    matrix JK0[`c', 1] = r(paf)
}
drop newid
mata: st_local("sejk0", strofreal(sqrt((rows(st_matrix("JK0")) - 1)/rows(st_matrix("JK0")) * ///
    sum((st_matrix("JK0") :- mean(st_matrix("JK0"))):^2)), "%21.15g"))
qui stcox z age, vce(cluster id) breslow
qui stpaf, eliminate(z) vce(jackknife)
local senw0 = r(se)
di as txt "no strata: IF SE = " %9.6f `seif0' "   vce(jackknife) SE = " %9.6f r(se) ///
    "   explicit loop SE = " %9.6f `sejk0'
assert reldif(`senw0', `sejk0') < 1e-5
assert inrange(`seif0'/`sejk0', 0.90, 1.05)
qui stpaf, eliminate(z) vce(jackknife) nowarm
assert reldif(r(se), `senw0') < 1e-9
di as txt "[ok] 10: vce(jackknife) without strata reproduces the explicit loop"

* ---- 11. cluster() option
* (a) vce(robust) on stset id() data: clusters come from the stset id
qui stcox z age, vce(robust) breslow
qui stpaf, eliminate(z)
assert reldif(r(paf), `paf0') < 1e-10 & reldif(r(se), `seif0') < 1e-10
assert "`r(clustvar)'" == "id"
* (b) cluster(id) gives the same
qui stpaf, eliminate(z) cluster(id)
assert reldif(r(paf), `paf0') < 1e-10 & reldif(r(se), `seif0') < 1e-10
assert r(N_clust) == `G'
* (c) coarser clusters (subjects nested): IF and jackknife run on 20 clusters
gen int grp = 1 + mod(id, 20)
qui stpaf, eliminate(z) cluster(grp)
assert r(N_clust) == 20 & reldif(r(paf), `paf0') < 1e-10
local se20 = r(se)
qui stpaf, eliminate(z) cluster(grp) vce(jackknife) nodots
assert r(N_clust) == 20
local sejk20 = r(se)
* (d) non-nested clusters (time-varying z): runs with a note, 2 clusters
stpaf, eliminate(z) cluster(z)
assert r(N_clust) == 2
* (e) missing cluster values are refused
qui replace grp = . if e(sample) & sum(e(sample)) == 1
capture stpaf, eliminate(z) cluster(grp)
assert _rc == 459
drop grp
* (f) data stset like ghoshlin output: [pw=weight_c], no id(), enter()
preserve
gen double weight_c = 1
stset t1b [pw = weight_c], failure(db) enter(t0b)
qui stcox z age, vce(robust) breslow
capture stpaf, eliminate(z)
assert _rc == 459
qui stpaf, eliminate(z) cluster(id)
assert reldif(r(paf), `paf0') < 1e-8 & reldif(r(se), `seif0') < 1e-8
qui stcox z age, vce(cluster id) breslow
qui stpaf, eliminate(z)
assert reldif(r(se), `seif0') < 1e-8
restore
di as txt "[ok] 11: cluster() (20 clusters: IF SE = " %6.4f `se20' ", jackknife SE = " %6.4f `sejk20' ///
    "); ghoshlin-type data need cluster() or vce(cluster)"

* ---- 7. weights (single-record Cox data, no id: each record is its own cluster)
clear
set seed 4711
set obs 300
gen double x = rnormal()
gen byte trt = runiform() > .5
gen double t = rexponential(1/exp(0.5*trt + 0.3*x))
gen byte d = runiform() < .8
gen byte fw = 1 + floor(3*runiform())
stset t [fw = fw], failure(d)
qui stcox trt x, breslow
qui stpaf, eliminate(trt)
local pf = r(paf)
local sf = r(se)
expand fw
stset t, failure(d)
qui stcox trt x, breslow
qui stpaf, eliminate(trt)
assert reldif(r(paf), `pf') < 1e-8
assert reldif(r(se), `sf') < 1e-6
di as txt "[ok] 7a: fweights reproduce the expanded data set"
gen byte one = 1
stset t [pw = one], failure(d)
qui stcox trt x, breslow
qui stpaf, eliminate(trt)
local pp = r(paf)
local sp = r(se)
stset t, failure(d)
qui stcox trt x, breslow
qui stpaf, eliminate(trt)
assert reldif(r(paf), `pp') < 1e-10 & reldif(r(se), `sp') < 1e-8
di as txt "[ok] 7b: unit pweights reproduce the unweighted result"

* ---- 12. ties
clear
set seed 815
set obs 300
gen long idv = _n
gen double x = rnormal()
gen byte trt = runiform() > .5
gen double t = rexponential(1/exp(0.5*trt + 0.3*x))
gen byte d = runiform() < .8
* (a) continuous times, efron: no tied failure times -> Newton, identical to breslow
stset t, id(idv) failure(d)
qui stcox trt x, efron vce(cluster idv)
qui stpaf, eliminate(trt) vce(jackknife)
assert "`r(jk_method)'" == "newton"
local sef = r(se)
qui stcox trt x, breslow vce(cluster idv)
qui stpaf, eliminate(trt) vce(jackknife)
assert reldif(r(se), `sef') < 1e-9
* (b) tied failure times, efron -> stcox refits (with a note)
gen double tr = ceil(t*4)/4
stset tr, id(idv) failure(d)
qui stcox trt x, efron vce(cluster idv)
stpaf, eliminate(trt) vce(jackknife) nodots
assert "`r(jk_method)'" == "refit"
* (c) tied failure times, breslow: Newton = refits
qui stcox trt x, breslow vce(cluster idv)
qui stpaf, eliminate(trt) vce(jackknife)
assert "`r(jk_method)'" == "newton"
local sen = r(se)
qui stpaf, eliminate(trt) vce(jackknife) refit nodots
assert reldif(r(se), `sen') < 1e-5
di as txt "[ok] 12: efron without ties -> Newton; efron with ties -> refit; breslow with ties: Newton = refit"

di as res _n "stpaf certification finished"
