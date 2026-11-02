* stpaf_ghoshlin_check.do -- compatibility checks of stpaf 0.4.0 with ghoshlin 1.3 (01nov2026)
* Needs stpaf.ado and ghoshlin.ado on the adopath. Uses ghoshlin_example.dta (961 subjects).
* Written without access to Stata: please run once and report the first failure, if any.
*
*  1. Ghosh-Lin data, vce(cluster ID): E = observed events, N_clust = subjects
*  2. Analytic case: PAF among exposed = 1 - exp(-b)
*  3. cluster() option: equals vce(cluster ID); guard after vce(robust); coarser clusters; missing values
*  4. IF vs vce(jackknife); Newton jackknife (warm, nowarm) vs stcox refits (set local refitcheck 0
*     to skip the slow refits)
*  5. Bootstrap with ghoshlin weights held fixed vs re-running ghoshlin in every replicate
*     (size of the omitted censoring-estimation term)
*  6. String ID gives identical results
*  7. fweights: default expansion (repid) = manually expanded data; noexpandfw same PAF
*  8. maxevents(1) (Fine-Gray data) runs; E = observed CV deaths
*  9. mi set data: not supported by stpaf (information only)

version 19
clear all
set more off
local B 200            // bootstrap replications for check 5 (500 for a sharper answer)
local refitcheck 0     // 1: check 4 also runs the stcox-refit jackknife (about 7-8 minutes here)
local seed 20261102

local stopt  failure(events == 1, 2, 3) exit(time .)
local glopt  failure("inlist(events,1,3)") competing("inlist(events,1,2)") keep(age sex diabetes)
local model  i.diabetes age i.sex

sysuse ghoshlin_example, clear
frame rename default orig

* ---------------------------------------------------------------- 1. basic run
frame copy orig work, replace
frame work {
    stset fup_days, id(ID) `stopt'
    ghoshlin, `glopt' frame(gl, replace)
}
frame gl {
    qui count if weight_c == 0
    di as txt "(rows with weight_c == 0: " r(N) ")"
    stcox `model', vce(cluster ID) breslow nolog
    stpaf, eliminate(diabetes)
    local paf1 = r(paf)
    local se1  = r(se)
    local E1   = r(E)
    local G1   = r(N_clust)
    qui summarize event [iw = weight_c] if e(sample), meanonly
    assert reldif(`E1', r(sum)) < 1e-8
    tempvar tag
    egen byte `tag' = tag(ID) if e(sample)
    qui count if `tag' == 1
    assert `G1' == r(N)
}
di as txt "[ok] 1: E = weighted observed events (" %8.2f `E1' "), clusters = subjects (" `G1' ")"
di as txt "        PAF = " %7.4f `paf1' "  SE(IF) = " %7.4f `se1'

* ---------------------------------------------------------------- 2. analytic case
frame gl {
    local bd = _b[1.diabetes]
    stpaf if diabetes == 1, eliminate(diabetes)
    assert reldif(r(paf), 1 - exp(-`bd')) < 1e-8
}
di as txt "[ok] 2: PAF among exposed = 1 - exp(-b)"

* ---------------------------------------------------------------- 3. cluster() option (stpaf 0.3.0)
frame gl {
    * (a) explicit cluster() equal to e(clustvar): identical results
    qui stcox `model', vce(cluster ID) breslow
    qui stpaf, eliminate(diabetes) cluster(ID)
    assert reldif(r(paf), `paf1') < 1e-10 & reldif(r(se), `se1') < 1e-10
    assert "`r(clustvar)'" == "ID"
    * (b) vce(robust): no e(clustvar), no stset id -> ghoshlin guard must stop stpaf
    qui stcox `model', vce(robust) breslow
    capture stpaf, eliminate(diabetes)
    assert _rc == 459
    * (c) vce(robust) + cluster(ID): same as fitting stcox with vce(cluster ID)
    qui stpaf, eliminate(diabetes) cluster(ID)
    assert reldif(r(paf), `paf1') < 1e-10 & reldif(r(se), `se1') < 1e-10
    assert r(N_clust) == `G1'
    * (d) jackknife follows cluster()
    qui stpaf, eliminate(diabetes) cluster(ID) vce(jackknife) nodots
    local sejk_c = r(se)
    * (e) coarser clusters (subjects nested): runs, fewer clusters
    gen int _grp = mod(ID, 50)
    qui stpaf, eliminate(diabetes) cluster(_grp)
    assert r(N_clust) == 50
    local se50 = r(se)
    * (f) missing cluster values are refused
    qui replace _grp = . if e(sample) & sum(e(sample)) == 1    // first estimation-sample row
    capture stpaf, eliminate(diabetes) cluster(_grp)
    assert _rc == 459
    drop _grp
}
di as txt "[ok] 3: cluster() = vce(cluster ID); guard after vce(robust); SE with 50 groups = " %7.4f `se50'

* ---------------------------------------------------------------- 4. IF vs jackknife
frame gl {
    qui stcox `model', vce(cluster ID) breslow
    timer clear
    timer on 1
    qui stpaf, eliminate(diabetes) vce(jackknife)
    timer off 1
    local sejk = r(se)
    local itw  = r(jk_iter)
    assert "`r(jk_method)'" == "newton"
    assert reldif(r(paf), `paf1') < 1e-10
    timer on 2
    qui stpaf, eliminate(diabetes) vce(jackknife) nowarm
    timer off 2
    local itc = r(jk_iter)
    assert reldif(r(se), `sejk') < 1e-9
    if `refitcheck' {
        timer on 3
        qui stpaf, eliminate(diabetes) vce(jackknife) refit nodots
        timer off 3
        local serf = r(se)
        assert "`r(jk_method)'" == "refit"
        assert reldif(`serf', `sejk') < 1e-5
    }
    * (stcox allows no weights with efron ties, so ghoshlin data always use the Newton jackknife;
    *  the efron fallback is tested in stpaf_cert.do, test 12)
}
qui timer list
di as txt "[info] 4: SE IF = " %7.4f `se1' "  SE jackknife = " %7.4f `sejk' ///
    "  ratio = " %5.3f `se1'/`sejk'
di as txt "         Newton warm start " %5.1f r(t1) " s (" %4.2f `itw' " it.), nowarm " %5.1f r(t2) ///
    " s (" %4.2f `itc' " it.)" _c
if `refitcheck' di as txt ", stcox refits " %5.1f r(t3) " s; SE refits = " %9.6f `serf' " vs Newton " %9.6f `sejk'
else di
assert inrange(`se1'/`sejk', 0.85, 1.15)
assert reldif(`sejk_c', `sejk') < 1e-8           // jackknife with cluster(ID) after vce(robust)
di as txt "[ok] 4"

* ---------------------------------------------------------------- 5. censoring-estimation term
tempname pf
tempfile bs
postfile `pf' int rep double(paf_fixed paf_full) using `bs'
set seed `seed'
forvalues r = 1/`B' {
    * (a) ghoshlin weights fixed: resample subjects of the transformed data
    frame copy gl bsf, replace
    frame bsf {
        qui bsample, cluster(ID) idcluster(newid)
        qui stset tstop [pw = weight_c], failure(event) enter(tstart)
        qui stcox `model', vce(cluster newid) breslow
        qui stpaf, eliminate(diabetes)
        local a = r(paf)
    }
    * (b) full: resample subjects of the original data and re-run ghoshlin
    frame copy orig bso, replace
    frame bso {
        qui bsample, cluster(ID) idcluster(newid)
        qui stset fup_days, id(newid) `stopt'
        qui ghoshlin, `glopt' frame(bsgl, replace)
    }
    frame bsgl {
        qui stcox `model', vce(cluster newid) breslow
        qui stpaf, eliminate(diabetes)
        local b = r(paf)
    }
    post `pf' (`r') (`a') (`b')
    if mod(`r', 50) == 0 di as txt "`r' " _c
}
postclose `pf'
frame create bsres
frame bsres {
    use `bs', clear
    qui summarize paf_fixed
    local sefx = r(sd)
    qui summarize paf_full
    local sefu = r(sd)
}
frame drop bsres bsf bso bsgl
di _n as txt "[info] 5: SE  IF = " %7.4f `se1' "   boot, weights fixed = " %7.4f `sefx' ///
    "   boot, ghoshlin re-run = " %7.4f `sefu'
di as txt "         IF/fixed = " %5.3f `se1'/`sefx' "   full/fixed = " %5.3f `sefu'/`sefx' ///
    "   (MC error of each ratio ~ " %4.2f 1/sqrt(2*`B') ")"
assert inrange(`se1'/`sefx', 0.85, 1.15)
di as txt "[ok] 5: IF agrees with the fixed-weight bootstrap"

* ---------------------------------------------------------------- 6. string ID
frame copy orig work, replace
frame work {
    qui tostring ID, replace format(%05.0f)
    qui stset fup_days, id(ID) `stopt'
    qui ghoshlin, `glopt' frame(gls, replace)
}
frame gls {
    qui stcox `model', vce(cluster ID) breslow
    qui stpaf, eliminate(diabetes)
    assert reldif(r(paf), `paf1') < 1e-8 & reldif(r(se), `se1') < 1e-8
    qui stpaf, eliminate(diabetes) vce(jackknife) nodots
    assert reldif(r(se), `sejk') < 1e-8
}
frame drop gls
di as txt "[ok] 6: string ID gives identical PAF, IF and jackknife SE"

* ---------------------------------------------------------------- 7. fweights
frame orig: gen byte fw = 1 + mod(ID, 3)          // constant within subject

* (a) default: ghoshlin expands into repid
frame copy orig work, replace
frame work {
    qui stset fup_days [fw = fw], id(ID) `stopt'
    qui ghoshlin, `glopt' frame(glfa, replace)
}
frame glfa {
    qui stcox `model', vce(cluster repid) breslow
    qui stpaf, eliminate(diabetes)
    local pa = r(paf)
    local sa = r(se)
}
* (b) manual expansion of the original data, unweighted
frame copy orig work, replace
frame work {
    gen long _o = _n
    qui expand fw
    bysort _o: gen int _copy = _n
    egen long newid = group(ID _copy)
    qui stset fup_days, id(newid) `stopt'
    qui ghoshlin, `glopt' frame(glfb, replace)
}
frame glfb {
    qui stcox `model', vce(cluster newid) breslow
    qui stpaf, eliminate(diabetes)
    local pb = r(paf)
    local sb = r(se)
}
* (c) noexpandfw
frame copy orig work, replace
frame work {
    qui stset fup_days [fw = fw], id(ID) `stopt'
    qui ghoshlin, `glopt' frame(glfc, replace) noexpandfw
}
frame glfc {
    qui stcox `model', vce(cluster ID) breslow
    qui stpaf, eliminate(diabetes)
    local pc = r(paf)
    local sc = r(se)
}
frame drop glfa glfb glfc
di as txt "[info] 7: PAF/SE expanded " %7.4f `pa' "/" %7.4f `sa' "   manual " %7.4f `pb' "/" %7.4f `sb' ///
    "   noexpandfw " %7.4f `pc' "/" %7.4f `sc'
assert reldif(`pa', `pb') < 1e-8 & reldif(`sa', `sb') < 1e-8
assert reldif(`pa', `pc') < 1e-8
di as txt "[ok] 7: expanded fweights = manual expansion; noexpandfw gives the same PAF"

* ---------------------------------------------------------------- 8. maxevents(1), Fine-Gray data
frame copy orig work, replace
frame work {
    qui stset fup_days, id(ID) failure(events == 1, 2) exit(time .)
    qui ghoshlin, failure("events==1") competing("events==2") keep(age sex diabetes) maxevents(1) frame(glfg, replace)
}
frame glfg {
    qui stcox `model', vce(cluster ID) breslow
    stpaf, eliminate(diabetes)
    local Efg = r(E)
    qui summarize event [iw = weight_c] if e(sample), meanonly
    assert reldif(`Efg', r(sum)) < 1e-8
}
frame drop glfg
di as txt "[ok] 8: maxevents(1) data; E = observed CV deaths"

* ---------------------------------------------------------------- 9. mi (information only)
di as txt "[info] 9: stpaf requires e(cmd2) == stcox; after -mi estimate: stcox- it stops by design."

di as res _n "stpaf / ghoshlin compatibility checks finished"
