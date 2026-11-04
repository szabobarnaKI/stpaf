* stpaf_cert_mi.do  --  certification of multiple-imputation support in stpaf 0.8.0 (test 16)
* Run with stpaf.ado on the adopath (e.g. -adopath ++ "path/to/folder"-), after stpaf_cert.do
* (which must still give 15/15 with 0.8.0).
*
* What is checked
*  16a. after mi estimate: stcox (influence function, 400 clusters): stpaf = explicit loop over
*       the imputations (mi extract, stcox, stpaf) combined by hand with Rubin's rules on
*       log(1-PAF), large-sample df (paf, se_logR, lb, ub, df to 1e-10); r(paf_mi) = loop
*  16b. the same with 20 clusters (default cluster jackknife): Barnard-Rubin df with complete
*       df 19; citrans(none): Rubin on the PAF scale
*  16c. the Rubin rules of 16a/16b (= stpaf) against mi estimate itself: large-sample df/FMI/total variance
*       (mi estimate: stcox) and small-sample df (mi estimate: regress, e(df_r))
*  16d. no between-imputation variability (model without the imputed variable): equals the
*       complete-data stpaf (IF: identical CI; jackknife: df = 19*20/22)
*  16e. mi estimate, post; mlong style; errors: data changed since mi estimate (498),
*       mi estimate of another command (301)

version 15.0
clear all
set more off
set type double
set seed 20261103

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
end

* per-imputation loop: stcox command `cmd', stpaf options `opt'; stores r(table) of each
* imputation stacked in matrix LOOP (M blocks)
capture program drop _cert_miloop
program define _cert_miloop
    syntax, cmd(string) opt(string) m(integer)
    capture matrix drop LOOP
    forvalues j = 1/`m' {
        preserve
        qui mi extract `j', clear
        qui `cmd'
        qui stpaf, `opt'
        matrix LOOP = nullmat(LOOP) \ r(table)
        restore
    }
end

* hand-coded Rubin's rules (independent of stpaf's Mata): matrix LOOP from _cert_miloop,
* scale lg = 1 (log) or 0 (PAF), complete df nuc (. = large sample); result matrix HAND
* (k x 7: paf, se, lb, ub, se_logR, df, fmi)
capture program drop _cert_rubin
program define _cert_rubin
    syntax, m(integer) lg(integer) nuc(string)
    mata: _cert_rubin_m(`m', `lg', `nuc')
end
mata:
// one column: estimate, total variance, df (missing = infinite), fmi (formulas of mi estimate)
real rowvector _cert_rub1(real colvector Q, real colvector U, real scalar nuc)
{
    real scalar M, qb, ub, B, T, r, nul, nuo, nu, fmi, i
    M = rows(Q)
    qb = mean(Q)
    ub = mean(U)
    B = 0
    for (i = 1; i <= M; i++) B = B + (Q[i] - qb)^2
    B = B/(M - 1)
    T = ub + (1 + 1/M)*B
    r = (1 + 1/M)*B/ub
    nul = (B > 0 ? (M - 1)*(1 + 1/r)^2 : .)
    if (nuc < .) {
        nuo = nuc*(nuc + 1)*(1 - (1 + 1/M)*B/T)/(nuc + 3)
        nu = (nul < . ? 1/(1/nul + 1/nuo) : nuo)
        fmi = 1 - ((nu + 1)/(nu + 3))/((nuc + 1)/(nuc + 3))*ub/T
    }
    else {
        nu = nul
        fmi = (nul < . ? (r + 2/(nul + 3))/(r + 1) : 0)
    }
    return((qb, T, nu, fmi))
}
void _cert_rubin_m(real scalar M, real scalar lg, real scalar nuc)
{
    real matrix L, H
    real scalar k, c, qb, T, nu, q
    real rowvector R
    real colvector Q, U, rows_c
    L = st_matrix("LOOP")
    k = rows(L)/M
    H = J(k, 7, .)
    for (c = 1; c <= k; c++) {
        rows_c = c :+ k*(0::M-1)
        if (lg) {
            Q = ln(1 :- L[rows_c, 3])
            U = L[rows_c, 9]:^2
        }
        else {
            Q = L[rows_c, 3]
            U = L[rows_c, 4]:^2
        }
        R = _cert_rub1(Q, U, nuc)
        qb = R[1]
        T = R[2]
        nu = R[3]
        q = (nu < . ? invttail(nu, .025) : invnormal(.975))
        if (lg) {
            H[c, .] = (1 - exp(qb), exp(qb)*sqrt(T), 1 - exp(qb + q*sqrt(T)),
                1 - exp(qb - q*sqrt(T)), sqrt(T), nu, R[4])
        }
        else {
            H[c, .] = (qb, sqrt(T), qb - q*sqrt(T), qb + q*sqrt(T), ., nu, R[4])
        }
    }
    st_matrix("HAND", H)
}
end

* ------------------------------------------------------------------ data
_cert_sim, n(400)
gen int grp = 1 + mod(id, 20)
* the data are stset by _cert_sim; clear the settings (mi would find stale st variables)
stset, clear
keep id t0b t1b db z age stratum grp
* age missing in about 30% of the subjects (all their records), then imputed
gen double u = runiform()
bysort id (t0b): replace u = u[1]
replace age = . if u < .3
drop u
qui count if missing(age)
assert r(N) > 0
mi set flong
mi register imputed age
mi register regular z stratum db t0b t1b grp id
mi impute regress age z i.stratum db, add(5) rseed(16)
mi stset t1b, id(id) failure(db) time0(t0b) exit(time .)
local cox "stcox z age, strata(stratum) vce(cluster id) breslow"

* ---- 16a. influence function (400 clusters), Rubin on log(1-PAF), large-sample df
qui mi estimate: `cox'
stpaf, eliminate(z) horizon(1 2.5)
assert r(M) == 5 & "`r(vce)'" == "if" & "`r(mi)'" == "mi"
assert r(df) < . & r(fmi) > 0 & r(fmi) < 1
matrix TM = r(table)
matrix PM = r(paf_mi)
assert colsof(TM) == 13 & rowsof(TM) == 2
_cert_miloop, cmd(`cox') opt(eliminate(z) horizon(1 2.5) vce(if)) m(5)
_cert_rubin, m(5) lg(1) nuc(.)
matrix A = TM[1..., 3..6], TM[1..., 9], TM[1..., 11], TM[1..., 13]
assert mreldif(A, HAND) < 1e-10
forvalues j = 1/5 {
    forvalues c = 1/2 {
        assert reldif(PM[`j', `c'], LOOP[(`j'-1)*2 + `c', 3]) < 1e-12
    }
}
di as txt "[ok] 16a: IF within imputations, Rubin on log(1-PAF): PAF " %6.4f TM[2, 3] ///
    " SE " %6.4f TM[2, 4] " df " %8.1f TM[2, 11] " FMI " %5.3f TM[2, 13]

* ---- 16b. 20 clusters: default jackknife within imputations, Barnard-Rubin df (complete df 19)
qui mi estimate: `cox'
stpaf, eliminate(z) horizon(1 2.5) cluster(grp)
assert "`r(vce)'" == "jackknife" & r(df_c) == 19 & r(N_clust) == 20
matrix TM = r(table)
_cert_miloop, cmd(`cox') opt(eliminate(z) horizon(1 2.5) cluster(grp) vce(jackknife)) m(5)
_cert_rubin, m(5) lg(1) nuc(19)
matrix A = TM[1..., 3..6], TM[1..., 9], TM[1..., 11], TM[1..., 13]
assert mreldif(A, HAND) < 1e-10
assert TM[1, 11] < 19 & TM[2, 11] < 19
local djk = TM[2, 11]
* citrans(none): Rubin on the PAF scale
qui mi estimate: `cox'
qui stpaf, eliminate(z) horizon(1 2.5) cluster(grp) citrans(none)
matrix TM = r(table)
_cert_rubin, m(5) lg(0) nuc(19)
matrix A = TM[1..., 3..6]
matrix H = HAND[1..., 1..4]
assert mreldif(A, H) < 1e-10
assert reldif(TM[2, 11], HAND[2, 6]) < 1e-10 & reldif(TM[2, 13], HAND[2, 7]) < 1e-10
di as txt "[ok] 16b: jackknife within imputations, Barnard-Rubin df " %6.2f `djk' " (complete df 19)"

* ---- 16c. Rubin's rules = mi estimate (df, total variance, FMI). The Mata functions of stpaf.ado
*          are private to the ado-file, so the hand-coded rules are checked here; 16a/16b tie
*          them to stpaf (estimate, se, CI, df, fmi to 1e-10)
* large sample: coefficient of age from stcox
qui mi estimate, vartable: `cox'
matrix DF = e(df_mi)
matrix VM = e(V_mi)
matrix FM = e(fmi_mi)
local ia 2    // columns z, age
local dfa = DF[1, `ia']
local va = VM[`ia', `ia']
local fa = FM[1, `ia']
matrix QU = J(5, 2, .)
forvalues j = 1/5 {
    preserve
    qui mi extract `j', clear
    qui `cox'
    matrix QU[`j', 1] = _b[age]
    matrix QU[`j', 2] = _se[age]^2
    restore
}
mata: R = _cert_rub1(st_matrix("QU")[., 1], st_matrix("QU")[., 2], .)
mata: st_local("rdf", strofreal(R[3], "%21.15g")); st_local("rt", strofreal(R[2], "%21.15g"))
mata: st_local("rfmi", strofreal(R[4], "%21.15g"))
assert reldif(`rdf', `dfa') < 1e-8 & reldif(`rt', `va') < 1e-8 & reldif(`rfmi', `fa') < 1e-8
* small sample: regress (stores e(df_r))
qui mi estimate, vartable: regress age z db
matrix DF = e(df_mi)
matrix VM = e(V_mi)
matrix FM = e(fmi_mi)
local dfz = DF[1, 1]
local vz = VM[1, 1]
local fz = FM[1, 1]
matrix QU = J(5, 2, .)
forvalues j = 1/5 {
    preserve
    qui mi extract `j', clear
    qui regress age z db
    matrix QU[`j', 1] = _b[z]
    matrix QU[`j', 2] = _se[z]^2
    local dfr = e(df_r)
    restore
}
mata: R = _cert_rub1(st_matrix("QU")[., 1], st_matrix("QU")[., 2], `dfr')
mata: st_local("rdf", strofreal(R[3], "%21.15g")); st_local("rt", strofreal(R[2], "%21.15g"))
mata: st_local("rfmi", strofreal(R[4], "%21.15g"))
di as txt "16c: regress z: mi estimate df " %9.3f `dfz' "  hand-coded " %9.3f `rdf' ///
    "   FMI " %6.4f `fz' " / " %6.4f `rfmi'
assert reldif(`rdf', `dfz') < 1e-8 & reldif(`rt', `vz') < 1e-8 & reldif(`rfmi', `fz') < 1e-8
di as txt "[ok] 16c: Rubin's rules and (Barnard-Rubin) df used by stpaf = mi estimate"

* ---- 16d. no between-imputation variability: equals the complete-data result
local cox0 "stcox z, strata(stratum) vce(cluster id) breslow"
qui mi estimate: `cox0'
qui stpaf, eliminate(z) vce(if)
matrix TM = r(table)
assert TM[1, 12] == 0 & TM[1, 13] == 0
local dfm = TM[1, 11]
preserve
qui mi extract 0, clear
qui `cox0'
qui stpaf, eliminate(z) vce(if)
matrix T0 = r(table)
restore
matrix A = TM[1, 3..6], TM[1, 9]
matrix B0 = T0[1, 3..6], T0[1, 9]
assert mreldif(A, B0) < 1e-10
assert `dfm' >= .
qui mi estimate: `cox0'
qui stpaf, eliminate(z) cluster(grp)
assert "`r(vce)'" == "jackknife"
assert reldif(r(df), 19*20/22) < 1e-10
local pj = r(paf)
local sj = r(se_logR)
preserve
qui mi extract 0, clear
qui `cox0'
qui stpaf, eliminate(z) cluster(grp)
assert reldif(r(paf), `pj') < 1e-10 & reldif(r(se_logR), `sj') < 1e-10
restore
di as txt "[ok] 16d: no imputed variable in the model: complete-data PAF and SE; jackknife df = 19*20/22"

* ---- 16e. post, mlong, errors
qui mi estimate: `cox'
qui stpaf, eliminate(z)
local p1 = r(paf)
local s1 = r(se)
qui mi estimate, post: `cox'
qui stpaf, eliminate(z)
assert r(M) == 5 & reldif(r(paf), `p1') < 1e-12 & reldif(r(se), `s1') < 1e-12
preserve
mi convert mlong, clear
qui mi estimate: `cox'
qui stpaf, eliminate(z)
assert reldif(r(paf), `p1') < 1e-10 & reldif(r(se), `s1') < 1e-10
restore
* data changed after mi estimate: the refits no longer reproduce e(b_mi)
qui mi estimate: `cox'
preserve
qui replace age = 2*age if _mi_m > 0
capture stpaf, eliminate(z)
assert _rc == 498
restore
* mi estimate of another command
qui mi estimate: regress age z db
capture stpaf, eliminate(z)
assert _rc == 301
di as txt "[ok] 16e: mi estimate, post and mlong give the same result; changed data -> 498; not stcox -> 301"

di as res _n "stpaf mi certification finished"
