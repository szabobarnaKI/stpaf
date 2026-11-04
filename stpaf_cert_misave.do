* stpaf_cert_misave.do  --  certification of stpaf 0.10.0: reuse of the stcox fits saved by
* mi estimate, saving() (test 18). Run with stpaf.ado on the adopath. Writes and erases
* f18*.ster files in the working directory.
*
* What is checked (saved = fits from the .ster file, refit = stcox refitted in each imputation)
*  18a. IF (400 clusters): saved path used in all 5 imputations (r(mi_nsaved) 5, r(mi_fits)
*       saved, r(mi_file)); r(table), r(paf_mi), r(se_mi) = mirefit = mi estimate without
*       saving() (refit fallback, nsaved 0); timings
*  18b. 20 clusters, default jackknife within imputations: saved = refit; timings
*  18c. mi estimate using f18a (stcox command taken from the file): = 18a;
*       using with nimputations(3): position mapping m, = mirefit = mi estimate, nimp(3)
*  18d. saving() with imputations(2 4) (file holds 2 fits): = mirefit; mifile() after an
*       mi estimate without saving() = 18a; mifile() of a missing file -> 601; saving file
*       erased -> note, refit, = 18a; mifile() with fits of another model -> refit, = refit
*  18e. mi estimate, post saving(); mlong style: = 18a
*  18f. data changed after mi estimate: saved fits fail the score check, refits do not
*       reproduce e(b_mi) -> 498 (also with mirefit)
*  18g. efron ties with tied failure times (Efron score check): saved = refit, IF and jackknife
*  18h. exactm ties with tied failure times: score check cannot pass -> all refitted
*       (r(mi_nsaved) 0), = mirefit

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

* largest relative difference between two matrices (equal missing cells count as equal;
* missing in only one of them or different shapes -> missing)
mata:
real scalar _c18_md(string scalar a, string scalar b)
{
    real matrix A, B, D
    A = st_matrix(a)
    B = st_matrix(b)
    if (rows(A) != rows(B) | cols(A) != cols(B)) return(.)
    if (any((A :>= .) :!= (B :>= .))) return(.)
    D = editmissing(abs(A - B) :/ (abs(B) :+ 1), 0)
    return(max(D))
}
end

* keep r(table), r(paf_mi), r(se_mi) of the last stpaf as T<s>, P<s>, S<s>
capture program drop _c18_keep
program define _c18_keep
    args s
    matrix T`s' = r(table)
    matrix P`s' = r(paf_mi)
    matrix S`s' = r(se_mi)
end

* compare the kept results s1 and s2 (all three matrices) to tolerance tol
capture program drop _c18_same
program define _c18_same
    args s1 s2 tol lab
    local mx 0
    foreach x in T P S {
        mata: st_numscalar("__c18", _c18_md("`x'`s1'", "`x'`s2'"))
        local d = scalar(__c18)
        if `d' >= . {
            di as err "  `lab': matrices `x'`s1' and `x'`s2' differ in shape or missing cells"
            exit 9
        }
        local mx = max(`mx', `d')
    }
    di as txt "  `lab': max rel. diff " %9.2e `mx'
    assert `mx' < `tol'
end

* data and imputations: as stpaf_cert_mi.do (age missing for 30% of the subjects)
capture program drop _c18_data
program define _c18_data
    syntax, n(integer) [round(real 0)]
    _cert_sim, n(`n')
    gen int grp = 1 + mod(id, 20)
    stset, clear
    if `round' > 0 {
        * times on a grid: many tied failure times
        replace t0b = round(t0b, `round')
        replace t1b = round(t1b, `round')
        drop if t1b <= t0b
    }
    keep id t0b t1b db z age stratum grp
    gen double u = runiform()
    bysort id (t0b): replace u = u[1]
    replace age = . if u < .3
    drop u
    mi set flong
    mi register imputed age
    mi register regular z stratum db t0b t1b grp id
    mi impute regress age z i.stratum db, add(5) rseed(16)
    mi stset t1b, id(id) failure(db) time0(t0b) exit(time .)
end

foreach f in f18a f18b f18d f18e f18f f18g f18h f18m f18x {
    capture erase `f'.ster
}

* ================================================================== data set 1
_c18_data, n(400)
local cox "stcox z age, strata(stratum) vce(cluster id) breslow"
local cox0 "stcox z, strata(stratum) vce(cluster id) breslow"
local o1 "eliminate(z) horizon(1 2.5)"

* ---- 18a. influence function (400 clusters)
qui mi estimate, saving(f18a, replace): `cox'
timer clear
timer on 1
stpaf, `o1'
timer off 1
assert r(M) == 5 & r(mi_nsaved) == 5 & "`r(mi_fits)'" == "saved" & "`r(vce)'" == "if"
assert "`r(mi_file)'" == "f18a.ster"
_c18_keep A
timer on 2
qui stpaf, `o1' mirefit
timer off 2
assert r(mi_nsaved) == 0 & "`r(mi_fits)'" == "refit" & "`r(mi_file)'" == "f18a.ster"
_c18_keep AR
_c18_same A AR 1e-8 "18a saved vs mirefit"
qui mi estimate: `cox'
qui stpaf, `o1'
assert r(mi_nsaved) == 0 & "`r(mi_fits)'" == "refit" & "`r(mi_file)'" == ""
_c18_keep AN
_c18_same A AN 1e-8 "18a saved vs no saving()"
qui timer list
di as txt "[ok] 18a: IF, saved fits used in 5/5 imputations = refits; time saved " ///
    %6.2f r(t1) " s, refit " %6.2f r(t2) " s"

* ---- 18b. 20 clusters: default jackknife within imputations
qui mi estimate, saving(f18b, replace): `cox'
timer clear
timer on 1
qui stpaf, `o1' cluster(grp)
timer off 1
assert r(mi_nsaved) == 5 & "`r(vce)'" == "jackknife" & r(df_c) == 19
local jm "`r(jk_method)'"
_c18_keep B
timer on 2
qui stpaf, `o1' cluster(grp) mirefit
timer off 2
assert r(mi_nsaved) == 0 & "`r(jk_method)'" == "`jm'"
_c18_keep BR
_c18_same B BR 1e-8 "18b saved vs mirefit"
qui timer list
di as txt "[ok] 18b: jackknife (`jm'), saved = refits; time saved " %6.2f r(t1) ///
    " s, refit " %6.2f r(t2) " s"

* ---- 18c. mi estimate using
qui mi estimate using f18a
qui stpaf, `o1'
assert r(M) == 5 & r(mi_nsaved) == 5
_c18_keep C
_c18_same A C 1e-10 "18c using = 18a"
qui mi estimate using f18a, nimputations(3)
qui stpaf, `o1'
assert r(M) == 3 & r(mi_nsaved) == 3 & "`r(m_mi)'" == "1 2 3"
_c18_keep C3
qui stpaf, `o1' mirefit
assert r(mi_nsaved) == 0
_c18_keep C3R
_c18_same C3 C3R 1e-8 "18c using nimp(3), saved vs mirefit"
qui mi estimate, nimputations(3): `cox'
qui stpaf, `o1'
_c18_keep C3N
_c18_same C3 C3N 1e-8 "18c using nimp(3) vs mi estimate, nimp(3)"
di as txt "[ok] 18c: mi estimate using (all and nimputations(3)) = refits"

* ---- 18d. subsets of imputations, mifile(), missing file, fits of another model
qui mi estimate, saving(f18d, replace) imputations(2 4): `cox'
* diagnostics: what the file holds (numbering of the fits of mi estimate's saving())
di as txt "  e(m_est_mi) = `e(m_est_mi)'"
estimates describe using f18d
forvalues k = 1/5 {
    capture quietly estimates use f18d, number(`k')
    local rck = _rc
    capture confirm matrix e(b)
    if `rck' di as txt "  estimates use f18d, number(`k'): rc = `rck'"
    else if _rc {
        di as txt "  estimates use f18d, number(`k'): no e(b); e(cmd) = `e(cmd)', e(cmd2) = `e(cmd2)'; ereturn list:"
        ereturn list
    }
    else di as txt "  estimates use f18d, number(`k'): b[age] = " %12.8f _b[age] ", N = " e(N)
}
* b[age] of imputations 2 and 4 refitted, to see which number holds which imputation
forvalues m = 2(2)4 {
    preserve
    qui mi extract `m', clear
    qui `cox'
    di as txt "  refit imputation `m': b[age] = " %12.8f _b[age]
    restore
}
* the diagnostics replaced e(): rerun mi estimate
qui mi estimate, saving(f18d, replace) imputations(2 4): `cox'
stpaf, `o1'
di as txt "  r(M) = " r(M) ", r(mi_nsaved) = " r(mi_nsaved) ", r(m_mi) = `r(m_mi)'"
assert r(M) == 2
assert "`r(m_mi)'" == "2 4"
assert r(mi_nsaved) == 2
_c18_keep D
qui stpaf, `o1' mirefit
_c18_keep DR
_c18_same D DR 1e-8 "18d imputations(2 4), saved vs mirefit"
* mi estimate without saving(), fits of 18a given in mifile()
qui mi estimate: `cox'
qui stpaf, `o1' mifile(f18a)
assert r(mi_nsaved) == 5 & "`r(mi_file)'" == "f18a.ster"
_c18_keep D2
_c18_same A D2 1e-10 "18d mifile() = 18a"
capture stpaf, `o1' mifile(f18nonexistent)
assert _rc == 601
* saving() file no longer there: note, refits
qui mi estimate, saving(f18x, replace): `cox'
erase f18x.ster
stpaf, `o1'
assert r(mi_nsaved) == 0 & "`r(mi_file)'" == ""
_c18_keep D3
_c18_same A D3 1e-8 "18d missing file -> refits = 18a"
* fits of another model (with age) given for a model without age: refits
qui mi estimate: `cox0'
qui stpaf, eliminate(z)
_c18_keep D4N
stpaf, eliminate(z) mifile(f18a)
assert r(mi_nsaved) == 0
_c18_keep D4
_c18_same D4 D4N 1e-10 "18d fits of another model -> refits"
di as txt "[ok] 18d: imputations(), mifile(), missing file (note), other model (note)"

* ---- 18e. post, mlong
qui mi estimate, post saving(f18e, replace): `cox'
qui stpaf, `o1'
assert r(mi_nsaved) == 5
_c18_keep E
_c18_same A E 1e-10 "18e post = 18a"
preserve
mi convert mlong, clear
qui mi estimate, saving(f18m, replace): `cox'
qui stpaf, `o1'
assert r(mi_nsaved) == 5
_c18_keep EM
_c18_same A EM 1e-8 "18e mlong = 18a"
restore
di as txt "[ok] 18e: post and mlong use the saved fits, = 18a"

* ---- 18f. data changed after mi estimate
qui mi estimate, saving(f18f, replace): `cox'
preserve
qui replace age = 2*age if _mi_m > 0
capture noisily stpaf, `o1'
assert _rc == 498
capture stpaf, `o1' mirefit
assert _rc == 498
restore
di as txt "[ok] 18f: changed data -> saved fits rejected (score), refits -> 498"

* ================================================================== data set 2: tied failure times
_c18_data, n(400) round(.25)
local coxe "stcox z age, strata(stratum) vce(cluster id) efron"

* ---- 18g. efron ties: the score check uses the Efron score
qui mi estimate, saving(f18g, replace): `coxe'
qui stpaf, `o1'
assert r(mi_nsaved) == 5 & "`r(vce)'" == "if"
_c18_keep G
qui stpaf, `o1' mirefit
_c18_keep GR
_c18_same G GR 1e-8 "18g efron, IF, saved vs mirefit"
qui stpaf, `o1' cluster(grp)
assert r(mi_nsaved) == 5 & "`r(vce)'" == "jackknife"
_c18_keep GJ
qui stpaf, `o1' cluster(grp) mirefit
_c18_keep GJR
_c18_same GJ GJR 1e-8 "18g efron, jackknife, saved vs mirefit"
di as txt "[ok] 18g: efron with tied failure times, saved = refits (IF and jackknife)"

* ---- 18h. exactm ties with tied failure times: saved fits cannot be checked -> refits
* (stcox refuses exactm with vce(cluster): clusters given to stpaf)
local coxm "stcox z age, strata(stratum) exactm"
qui mi estimate, saving(f18h, replace): `coxm'
stpaf, `o1' cluster(id) vce(if)
assert r(mi_nsaved) == 0 & "`r(mi_fits)'" == "refit" & "`r(mi_file)'" == "f18h.ster"
_c18_keep H
qui stpaf, `o1' cluster(id) vce(if) mirefit
_c18_keep HR
_c18_same H HR 1e-10 "18h exactm, fallback vs mirefit"
di as txt "[ok] 18h: exactm with tied failure times -> refitted in all imputations (note)"

foreach f in f18a f18b f18d f18e f18f f18g f18h f18m {
    capture erase `f'.ster
}
di as res _n "stpaf saved-fits certification (test 18) finished"
