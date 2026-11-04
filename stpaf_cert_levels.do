* stpaf_cert_levels.do  --  certification of varlevels() (level-specific PAFs) in stpaf 0.11.0 (test 19)
* Run with stpaf.ado on the adopath (e.g. -adopath ++ "path/to/folder"-), or from run_stpaf_0110.do.
*
* Reference used throughout: moving level k of x to the reference r is the same as
* eliminate(dk) in the equivalent model with indicators dk = (x == k) and their interactions
* (same partial likelihood); the indicator path uses the scenario code validated in tests 1-18.
*
* What is checked
*  19a. (0) PAF and E_cf of varlevels(x), set(x=0) and the indicator model = computed by hand
*       (Mata, columns built directly from x); then
*       varlevels(x) after stcox i.x##c.age = eliminate(d1 d2 d3) each after the indicator model,
*       whole r(table) (2 horizons x 4 scenarios): influence function (400 clusters), default
*       Newton jackknife (20 clusters), stcox refits; efron ties with tied failure times (Efron IF
*       and Efron Newton jackknife). Data: recurrent events, 2 strata, delayed entry, x time-varying
*       within subject (levels 0-3), interaction with age.
*  19b. additivity: sum of the level PAFs = joint PAF and sum of (E - E_cf) = joint (E - E_cf) at
*       each horizon (1e-12), with an if subset; r(sumdiff); joint row = set(x=0) (IF and jackknife)
*  19c. reference: ib2.x model (default reference 2) = i.x model with varlevels(x=2) = indicator
*       model with reference 2; a level without records in the if subset gives PAF 0, SE 0
*       (IF and jackknife); ibn.x without a reference -> 198
*  19d. errors: x continuous in the model, variable not in the model as a factor, combination with
*       set()/each, reference not a level, non-integer reference, unknown variable
*  19e. after mi estimate (age imputed): varlevels(x) (saved fits) = eliminate(d1 d2 d3) each after
*       mi estimate of the indicator model (refits), IF and jackknife; per-imputation level PAFs add
*       up exactly (r(paf_mi)); citrans(none): combined level PAFs add up (r(sumdiff) < 1e-12);
*       citrans(log): r(sumdiff) = hand-computed difference

version 15.0
clear all
set more off
set type double
set seed 20261104

* data: n subjects, recurrent events, 2 strata, 30% delayed entry; baseline level xb of x in
* {0,1,2}, one level up after the switch time sw (x time-varying, levels 0-3); round: event times
* rounded up to multiples of 0.25 (tied failure times)
capture program drop _clv_sim
program define _clv_sim
    syntax, n(integer) [ROUND]
    clear
    set obs `n'
    gen long id = _n
    gen double age = rnormal()
    gen byte stratum = 1 + (runiform() > .5)
    gen byte xb = floor(3*runiform())
    gen double C = 2 + 3*runiform()
    gen double sw = rexponential(2)
    gen double fr = rgamma(2, .5)
    gen double e0 = cond(runiform() < .3, .5*runiform(), 0)
    gen double len0 = min(sw, C)
    gen double len1 = max(C - sw, 0)
    gen double r0 = (0.6 + 0.4*(stratum == 2))*fr*exp(0.3*age + 0.35*xb + 0.1*xb*age)
    gen double r1 = (0.6 + 0.4*(stratum == 2))*fr*exp(0.3*age + 0.35*(xb + 1) + 0.1*(xb + 1)*age)
    gen long n0 = rpoisson(r0*len0)
    gen long n1 = rpoisson(r1*len1)
    expand n0 + n1 + 1
    bysort id: gen long k = _n
    gen double tev = .
    replace tev = runiform()*len0 if k <= n0
    replace tev = sw + runiform()*len1 if k > n0 & k <= n0 + n1
    replace tev = C if k == n0 + n1 + 1
    gen byte ev = k <= n0 + n1
    if "`round'" != "" replace tev = ceil(4*tev)/4
    gsort id tev -ev
    by id: gen double tstart = cond(_n == 1, 0, tev[_n-1])
    drop if tev <= tstart
    stset tev, id(id) failure(ev) time0(tstart) enter(time e0) exit(time .)
    stsplit blk, after(time = sw) at(0)
    gen byte x = xb + (_t0 >= sw)
    forvalues j = 0/3 {
        gen byte d`j' = (x == `j')
    }
    gen int grp = 1 + mod(id, 20)
end

* largest relative difference between two matrices (|a-b|/(|b|+1)); both missing counts as equal
capture mata: mata drop _clv_mrd()
mata:
real scalar _clv_mrd(string scalar an, string scalar bn)
{
    real matrix A, B, ma, mb
    A = st_matrix(an)
    B = st_matrix(bn)
    if (rows(A) != rows(B) | cols(A) != cols(B)) return(1e300)
    ma = (A :>= .)
    mb = (B :>= .)
    if (any(ma :!= mb)) return(1e300)
    A = editmissing(A, 0)
    B = editmissing(B, 0)
    return(max(abs(A - B) :/ (abs(B) :+ 1)))
}
end
capture program drop _clv_cmp
program define _clv_cmp
    args a b tol what
    mata: st_numscalar("CLVD", _clv_mrd("`a'", "`b'"))
    di as txt "   `what': max rel. difference " %9.2e scalar(CLVD)
    assert scalar(CLVD) < `tol'
end

* additivity check of an r(table) T with nh horizons: sum of level PAFs and of E - E_cf
capture program drop _clv_add
program define _clv_add
    args T nh what
    local ns = rowsof(`T')/`nh'
    local dmax 0
    forvalues h = 1/`nh' {
        scalar SP = 0
        scalar SE = 0
        forvalues s = 2/`ns' {
            local r = (`s' - 1)*`nh' + `h'
            scalar SP = SP + `T'[`r', 3]
            scalar SE = SE + (`T'[`r', 7] - `T'[`r', 8])
            assert reldif(`T'[`r', 7], `T'[`h', 7]) < 1e-12
        }
        local d1 = abs(SP - `T'[`h', 3])
        local d2 = reldif(SE, `T'[`h', 7] - `T'[`h', 8])
        local dmax = max(`dmax', `d1', `d2')
    }
    di as txt "   `what': additivity, max difference " %9.2e `dmax'
    assert `dmax' < 1e-12
end

* PAF by hand (Breslow, strata, (t0, t1] records, no tied failure times) for the i.x##c.age model
* in memory: scenario 1 all levels of x to 0, scenario 1+k level k only; horizons 1 and 2.5;
* indicator and interaction columns built directly from x, independently of stpaf and fvrevar;
* result matrix HAND (8 x 4: scenario, horizon, paf, E_cf)
capture program drop _clv_hand
program define _clv_hand
    tempvar es
    gen byte `es' = e(sample)
    matrix HB = (_b[1.x], _b[2.x], _b[3.x], _b[age], _b[1.x#c.age], _b[2.x#c.age], _b[3.x#c.age])
    forvalues k = 1/3 {
        tempvar h`k' g`k'
        gen double `h`k'' = (x == `k')
        gen double `g`k'' = (x == `k')*age
    }
    mata: _clv_hand_m("`h1' `h2' `h3' age `g1' `g2' `g3'", "`es'")
end
capture mata: mata drop _clv_hand_m()
mata:
void _clv_hand_m(string scalar xl, string scalar es)
{
    real matrix X, Xc, AC, H
    real colvector xv, t0, t1, dd, st, a, m, ev, risk
    real rowvector b, E0, tau
    real scalar s, k, ti, i, r, tt, S0, E
    b = st_matrix("HB")
    X = st_data(., tokens(xl), es)
    xv = st_data(., "x", es)
    t0 = st_data(., "_t0", es)
    t1 = st_data(., "_t", es)
    dd = st_data(., "_d", es)
    st = st_data(., "stratum", es)
    a = exp(X*b')
    AC = J(rows(X), 4, .)
    for (s = 1; s <= 4; s++) {
        Xc = X
        for (k = 1; k <= 3; k++) {
            if (s == 1 | s == k + 1) {
                m = (xv :== k)
                Xc[., k] = Xc[., k] :* (1 :- m)
                Xc[., k + 4] = Xc[., k + 4] :* (1 :- m)
            }
        }
        AC[., s] = exp(Xc*b')
    }
    ev = selectindex(dd :== 1)
    H = J(8, 4, .)
    tau = (1, 2.5)
    for (ti = 1; ti <= 2; ti++) {
        E = 0
        E0 = J(1, 4, 0)
        for (i = 1; i <= rows(ev); i++) {
            r = ev[i]
            if (t1[r] > tau[ti]) continue
            tt = t1[r]
            risk = selectindex((st :== st[r]) :& (t0 :< tt) :& (t1 :>= tt))
            S0 = sum(a[risk])
            E = E + 1
            E0 = E0 + colsum(AC[risk, .])/S0
        }
        for (s = 1; s <= 4; s++) H[(s - 1)*2 + ti, .] = (s, tau[ti], 1 - E0[s]/E, E0[s])
    }
    st_matrix("HAND", H)
}
end

* =========================================================================== 19a
_clv_sim, n(400)
local coxA "stcox i.x##c.age, strata(stratum) vce(cluster id) breslow"
local coxB "stcox d1 d2 d3 age c.d1#c.age c.d2#c.age c.d3#c.age, strata(stratum) vce(cluster id) breslow"
tab x if _st == 1
qui `coxB'
matrix BB = e(b)
`coxA'
matrix BA = e(b)
di as txt "19a: coefficient of 1.x: i.x model " %10.7f _b[1.x] ", indicator model " %10.7f BB[1, 1]
assert reldif(_b[1.x], BB[1, 1]) < 1e-6

* (0) point estimates by hand (independent of stpaf and fvrevar): varlevels(x), set(x=0) and the
*     indicator model with eliminate() (whose c.dk#c.age columns were stale before 0.11.0)
_clv_hand
matrix HP = HAND[1..., 3..4]
qui stpaf, varlevels(x) horizon(1 2.5) vce(if)
matrix T = r(table)
matrix TP = T[1..., 3], T[1..., 8]
_clv_cmp TP HP 1e-10 "19a(0) varlevels(x) PAF, E_cf vs hand"
qui stpaf, set(x=0) horizon(1 2.5) vce(if)
matrix T = r(table)
matrix TP = T[1..., 3], T[1..., 8]
matrix HJ = HP[1..2, 1...]
_clv_cmp TP HJ 1e-10 "19a(0) set(x=0) vs hand"
qui `coxB'
qui stpaf, eliminate(d1 d2 d3) each horizon(1 2.5) vce(if)
matrix T = r(table)
matrix TP = T[1..., 3], T[1..., 8]
_clv_cmp TP HP 1e-6 "19a(0) indicator model eliminate() each vs hand"
qui `coxA'

* (i) influence function, 400 clusters
stpaf, varlevels(x) horizon(1 2.5) vce(if)
assert "`r(levelvar)'" == "x" & r(reflevel) == 0 & "`r(levels)'" == "1 2 3"
assert "`r(scen1)'" == "x=0" & "`r(scen2)'" == "1.x" & "`r(scen3)'" == "2.x" & "`r(scen4)'" == "3.x"
assert "`r(vce)'" == "if" & r(N_clust) >= 390
matrix TA = r(table)
assert rowsof(TA) == 8 & colsof(TA) == 10
assert reldif(r(paf), TA[2, 3]) < 1e-15
qui `coxB'
qui stpaf, eliminate(d1 d2 d3) each horizon(1 2.5) vce(if)
matrix TB = r(table)
_clv_cmp TA TB 1e-6 "19a(i) IF, varlevels(x) vs eliminate(d1 d2 d3) each"

* (ii) default Newton jackknife, 20 clusters
qui `coxA'
stpaf, varlevels(x) horizon(1 2.5) cluster(grp)
assert "`r(vce)'" == "jackknife" & "`r(jk_method)'" == "newton" & r(df) == 19
matrix TA = r(table)
qui `coxB'
qui stpaf, eliminate(d1 d2 d3) each horizon(1 2.5) cluster(grp)
matrix TB = r(table)
_clv_cmp TA TB 1e-6 "19a(ii) Newton jackknife"
_clv_add TA 2 "19a(ii)"

* (iii) stcox refits
qui `coxA'
qui stpaf, varlevels(x) horizon(1 2.5) cluster(grp) refit nodots
assert "`r(jk_method)'" == "refit"
matrix TR = r(table)
_clv_cmp TR TA 1e-5 "19a(iii) refit vs Newton (i.x model)"
qui `coxB'
qui stpaf, eliminate(d1 d2 d3) each horizon(1 2.5) cluster(grp) refit nodots
matrix TB = r(table)
_clv_cmp TR TB 1e-5 "19a(iii) refit, i.x vs indicator model"

* (iv) efron ties with tied failure times
preserve
_clv_sim, n(300) round
local coxAe "stcox i.x##c.age, strata(stratum) vce(cluster id) efron"
local coxBe "stcox d1 d2 d3 age c.d1#c.age c.d2#c.age c.d3#c.age, strata(stratum) vce(cluster id) efron"
qui `coxAe'
qui stpaf, varlevels(x) horizon(1 2.5) vce(if)
assert "`r(if_ties)'" == "efron"
matrix TA = r(table)
qui `coxBe'
qui stpaf, eliminate(d1 d2 d3) each horizon(1 2.5) vce(if)
matrix TB = r(table)
_clv_cmp TA TB 1e-6 "19a(iv) efron IF"
qui `coxAe'
qui stpaf, varlevels(x) horizon(1 2.5) cluster(grp)
assert "`r(vce)'" == "jackknife" & "`r(jk_method)'" == "newton"
matrix TA = r(table)
qui `coxBe'
qui stpaf, eliminate(d1 d2 d3) each horizon(1 2.5) cluster(grp)
matrix TB = r(table)
_clv_cmp TA TB 1e-6 "19a(iv) efron Newton jackknife"
_clv_add TA 2 "19a(iv)"
restore
di as txt "[ok] 19a: varlevels() = indicator model with eliminate() each (IF, jackknife, refit, efron)"

* =========================================================================== 19b
qui `coxA'
qui stpaf if age > 0, varlevels(x) horizon(1 2.5) vce(if)
matrix T = r(table)
local sd = r(sumdiff)
_clv_add T 2 "19b IF, if age > 0"
assert `sd' < 1e-12
matrix TJ = T[1..2, 1...]
qui stpaf if age > 0, set(x=0) horizon(1 2.5) vce(if)
matrix S = r(table)
_clv_cmp TJ S 1e-12 "19b joint row = set(x=0), IF"
qui stpaf if age > 0, varlevels(x) horizon(1 2.5) cluster(grp)
matrix T = r(table)
_clv_add T 2 "19b jackknife, if age > 0"
matrix TJ = T[1..2, 1...]
qui stpaf if age > 0, set(x=0) horizon(1 2.5) cluster(grp)
matrix S = r(table)
_clv_cmp TJ S 1e-12 "19b joint row = set(x=0), jackknife"
di as txt "[ok] 19b: level PAFs add up to the joint PAF; joint row = set(x=0)"

* =========================================================================== 19c
* reference 2: ib2.x model (default), i.x model with varlevels(x=2), indicator model
qui stcox ib2.x##c.age, strata(stratum) vce(cluster id) breslow
qui stpaf, varlevels(x) horizon(1 2.5) vce(if)
assert r(reflevel) == 2 & "`r(levels)'" == "0 1 3" & "`r(scen1)'" == "x=2" & "`r(scen2)'" == "0.x"
matrix T2 = r(table)
qui `coxA'
qui stpaf, varlevels(x=2) horizon(1 2.5) vce(if)
assert r(reflevel) == 2 & "`r(levels)'" == "0 1 3"
matrix TA = r(table)
_clv_cmp TA T2 1e-6 "19c ib2.x default vs i.x with varlevels(x=2)"
qui stcox d0 d1 d3 age c.d0#c.age c.d1#c.age c.d3#c.age, strata(stratum) vce(cluster id) breslow
qui stpaf, eliminate(d0 d1 d3) each horizon(1 2.5) vce(if)
matrix TB = r(table)
_clv_cmp TA TB 1e-6 "19c varlevels(x=2) vs indicator model, reference 2"
_clv_add TA 2 "19c reference 2"
* level 3 without records in the if subset: PAF 0, SE 0
qui `coxA'
qui stpaf if x != 3, varlevels(x) horizon(1 2.5) vce(if)
matrix T = r(table)
forvalues r = 7/8 {
    assert abs(T[`r', 3]) < 1e-12 & abs(T[`r', 4]) < 1e-10 & abs(T[`r', 9]) < 1e-10
}
_clv_add T 2 "19c level without records, IF"
qui stpaf if x != 3, varlevels(x) horizon(1 2.5) cluster(grp)
matrix T = r(table)
forvalues r = 7/8 {
    assert abs(T[`r', 3]) < 1e-12 & abs(T[`r', 4]) < 1e-10 & abs(T[`r', 9]) < 1e-10
}
* ibn.x: no base level -> the reference must be given
capture qui stcox ibn.x age, strata(stratum) vce(cluster id) breslow
if !_rc {
    local cn : colnames e(b)
    di as txt "   ibn.x model, columns: `cn'"
    capture stpaf, varlevels(x)
    assert _rc == 198
}
else di as txt "   (stcox ibn.x not fitted, rc = " _rc "; subtest skipped)"
di as txt "[ok] 19c: reference level from e(b) or given; level without records -> PAF 0, SE 0"

* =========================================================================== 19d
qui `coxA'
capture stpaf, varlevels(x) set(age=0)
assert _rc == 198
capture stpaf, varlevels(x) each
assert _rc == 198
capture stpaf, varlevels(x=7)
assert _rc == 198
capture stpaf, varlevels(x=1.5)
assert _rc == 198
capture stpaf, varlevels(age)
assert _rc == 198
capture stpaf, varlevels(nosuchvar)
assert _rc == 111
qui stcox x age, strata(stratum) vce(cluster id) breslow
capture stpaf, varlevels(x)
assert _rc == 198
* blanks around = are allowed
qui `coxA'
qui stpaf, varlevels(x = 0) vce(if)
assert r(reflevel) == 0
di as txt "[ok] 19d: errors"

* =========================================================================== 19e
_clv_sim, n(400)
keep if _st == 1
gen double t0b = _t0
gen double t1b = _t
gen byte db = _d
stset, clear
keep id t0b t1b db x d0 d1 d2 d3 age stratum grp
gen double u = runiform()
bysort id (t0b): replace u = u[1]
replace age = . if u < .3
drop u
mi set flong
mi register imputed age
mi register regular x d0 d1 d2 d3 stratum db t0b t1b grp id
mi impute regress age i.x i.stratum db, add(5) rseed(17)
mi stset t1b, id(id) failure(db) time0(t0b) exit(time .)
local coxA "stcox i.x##c.age, strata(stratum) vce(cluster id) breslow"
local coxB "stcox d1 d2 d3 age c.d1#c.age c.d2#c.age c.d3#c.age, strata(stratum) vce(cluster id) breslow"
capture erase f19e.ster

* (i) influence function: saved fits (i.x model) vs refits (indicator model)
qui mi estimate, saving(f19e, replace): `coxA'
stpaf, varlevels(x) horizon(1 2.5) vce(if)
assert r(M) == 5 & r(mi_nsaved) == 5 & "`r(levels)'" == "1 2 3" & r(reflevel) == 0
matrix TA = r(table)
matrix PA = r(paf_mi)
scalar SDL = r(sumdiff)
local sdl = scalar(SDL)
assert colsof(TA) == 13 & rowsof(TA) == 8
qui mi estimate: `coxB'
qui stpaf, eliminate(d1 d2 d3) each horizon(1 2.5) vce(if)
matrix TB = r(table)
_clv_cmp TA TB 1e-6 "19e(i) mi, IF"
* per-imputation level PAFs add up exactly (columns of r(paf_mi) = rows of r(table))
scalar DM = 0
forvalues j = 1/5 {
    forvalues h = 1/2 {
        scalar SP = PA[`j', 2 + `h'] + PA[`j', 4 + `h'] + PA[`j', 6 + `h']
        scalar DM = max(DM, abs(SP - PA[`j', `h']))
    }
}
di as txt "   19e: per-imputation additivity, max difference " %9.2e scalar(DM)
assert scalar(DM) < 1e-12
* citrans(log): r(sumdiff) = difference of the combined estimates (approximate additivity)
scalar DL = 0
forvalues h = 1/2 {
    scalar SP = TA[2 + `h', 3] + TA[4 + `h', 3] + TA[6 + `h', 3]
    scalar DL = max(DL, abs(SP - TA[`h', 3]))
}
di as txt "   19e: citrans(log): sum of level PAFs - joint PAF up to " %9.2e scalar(DL) ///
    " (r(sumdiff) " %9.2e `sdl' ")"
assert reldif(`sdl', scalar(DL)) < 1e-6
* citrans(none): exact
qui mi estimate using f19e
qui stpaf, varlevels(x) horizon(1 2.5) vce(if) citrans(none)
assert r(mi_nsaved) == 5
assert r(sumdiff) < 1e-12
matrix T = r(table)
_clv_add T 2 "19e(i) citrans(none)"

* (ii) jackknife (20 clusters), Barnard-Rubin df
qui mi estimate, saving(f19e, replace): `coxA'
qui stpaf, varlevels(x) horizon(1 2.5) cluster(grp)
assert "`r(vce)'" == "jackknife" & r(df_c) == 19 & r(mi_nsaved) == 5
matrix TA = r(table)
qui mi estimate: `coxB'
qui stpaf, eliminate(d1 d2 d3) each horizon(1 2.5) cluster(grp)
matrix TB = r(table)
_clv_cmp TA TB 1e-6 "19e(ii) mi, jackknife"
capture erase f19e.ster
di as txt "[ok] 19e: varlevels() after mi estimate = indicator model; additivity per imputation and with citrans(none)"

di as res _n "stpaf varlevels() certification finished (19a-e)"
