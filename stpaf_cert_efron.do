* stpaf_cert_efron.do -- certification of the Efron-exact influence function (stpaf 0.9.0)
* Test 17: r(V_b) (cluster-robust VCE of b used by the influence function) must equal
* stcox's own vce(cluster) e(V), which Stata builds from its Efron score residuals and
* e(V_modelbased); without tied failure times the Breslow path is kept.
* Self-contained: own data generator (recurrent events, delayed entry, 2 strata, times
* rounded to a grid -> many tied failure times).
* Expected: all lines "PASS"; the run stops at the first failure.

version 15
clear all
set more off
discard
which stpaf
local t0 = clock("`c(current_time)'", "hms")

capture program drop _cef_sim
program define _cef_sim
    syntax , n(integer) seed(integer) [step(real 0.25)]
    clear
    set seed `seed'
    quietly {
        set obs `n'
        gen long id = _n
        gen double age = rnormal()
        gen byte z = runiform() < .4
        gen byte str = 1 + (runiform() < .5)
        gen double fr = rgamma(2, .5)
        gen double C = runiform(2, 5)
        gen double e0 = cond(runiform() < .3, runiform(), 0)
        gen double rate = .4*fr*exp(.7*z + .3*age)*cond(str == 2, 1.3, 1)
        expand 60
        bysort id: gen int j = _n
        by id: gen double t1 = e0 + sum(rexponential(1/rate))
        by id: gen double t0 = cond(_n == 1, e0, t1[_n-1])
        drop if t0 >= C
        gen byte d = t1 < C
        replace t1 = C if !d
        if `step' > 0 {
            replace t0 = round(t0, `step')
            replace t1 = round(t1, `step')
            drop if t1 <= t0
        }
        gen int grp = mod(id, 20) + 1
        stset t1, id(id) failure(d) time0(t0) exit(time .)
    }
end

capture program drop _cef_chk
program define _cef_chk
    args label val tol
    if `val' <= `tol' di as res "PASS" as txt "  `label' (" %9.2e `val' ")"
    else {
        di as err "FAIL  `label' (" %9.2e `val' " > `tol')"
        exit 9
    }
end

* ----------------------------------------------------------------------------- data
_cef_sim, n(400) seed(1701)
tempvar nf
quietly bysort str _t: egen `nf' = total(_d)
quietly summarize `nf' if _d, meanonly
di as txt "records " _N ", failures " r(N) ", largest number of tied failures " r(max)
assert r(max) > 1

* ----------------------------------------------------------------------------- 17a
di _n as txt "17a: efron ties, vce(cluster id): r(V_b) = e(V), r(if_ties) = efron"
stcox z age, efron strata(str) vce(cluster id)
tempname V1 TA VB
matrix `V1' = e(V)
stpaf, eliminate(z) horizon(2 5) vce(if)
assert "`r(if_ties)'" == "efron"
matrix `VB' = r(V_b)
matrix `TA' = r(table)
local m = mreldif(`VB', `V1')
_cef_chk "r(V_b) = stcox efron vce(cluster id) e(V)" `m' 1e-8
* the Efron-score check of e(b) (note about efron/exactp ties) must not fire: score = 0
di as txt "  SE(z) stcox " %9.6f sqrt(`V1'[1,1]) "  stpaf " %9.6f sqrt(`VB'[1,1])

* ----------------------------------------------------------------------------- 17b
di _n as txt "17b: same model without vce(), cluster(id) given to stpaf: same r(V_b) and r(table)"
stcox z age, efron strata(str)
stpaf, eliminate(z) horizon(2 5) vce(if) cluster(id)
assert "`r(if_ties)'" == "efron"
tempname VBb TAb
matrix `VBb' = r(V_b)
matrix `TAb' = r(table)
local m = mreldif(`VBb', `V1')
_cef_chk "r(V_b) = 17a" `m' 1e-8
local m = mreldif(`TAb', `TA')
_cef_chk "r(table) = 17a" `m' 1e-10

* ----------------------------------------------------------------------------- 17c
di _n as txt "17c: breslow ties: Breslow path, r(V_b) = e(V)"
stcox z age, breslow strata(str) vce(cluster id)
tempname V2
matrix `V2' = e(V)
stpaf, eliminate(z) horizon(2 5) vce(if)
assert "`r(if_ties)'" == "breslow"
tempname VBc
matrix `VBc' = r(V_b)
local m = mreldif(`VBc', `V2')
_cef_chk "r(V_b) = stcox breslow vce(cluster id) e(V)" `m' 1e-8

* ----------------------------------------------------------------------------- 17d
di _n as txt "17d: efron, clusters above the subject (20 groups): r(V_b) = stcox vce(cluster grp)"
stcox z age, efron strata(str) vce(cluster grp)
tempname V3
matrix `V3' = e(V)
stpaf, eliminate(z) vce(if)
assert "`r(if_ties)'" == "efron"
tempname VBd
matrix `VBd' = r(V_b)
local m = mreldif(`VBd', `V3')
_cef_chk "r(V_b) = stcox efron vce(cluster grp) e(V)" `m' 1e-8

* ----------------------------------------------------------------------------- 17e
di _n as txt "17e: efron, 400 clusters: influence function vs Newton (Efron) jackknife"
stcox z age, efron strata(str) vce(cluster id)
stpaf, eliminate(z) horizon(2 5) vce(jackknife)
tempname TJ
matrix `TJ' = r(table)
forvalues i = 1/2 {
    local ratio = `TA'[`i', 9]/`TJ'[`i', 9]
    di as txt "  horizon " %5.2f `TA'[`i', 2] ": se_logR IF " %8.5f `TA'[`i', 9] ///
        "  jackknife " %8.5f `TJ'[`i', 9] "  ratio " %6.3f `ratio'
    local m = abs(`ratio' - 1)
    _cef_chk "IF/jackknife ratio within 10%" `m' 0.10
}
tempname PA PJ
matrix `PA' = `TA'[1..., 3]
matrix `PJ' = `TJ'[1..., 3]
local m = mreldif(`PA', `PJ')
_cef_chk "point estimates IF path = jackknife path" `m' 1e-12

* ----------------------------------------------------------------------------- 17f
di _n as txt "17f: no tied failure times: efron model uses the Breslow path, r(V_b) = e(V)"
_cef_sim, n(300) seed(1702) step(0)
tempvar nf2
quietly bysort str _t: egen `nf2' = total(_d)
quietly summarize `nf2' if _d, meanonly
assert r(max) == 1
stcox z age, efron strata(str) vce(cluster id)
tempname V4
matrix `V4' = e(V)
stpaf, eliminate(z) vce(if)
assert "`r(if_ties)'" == "breslow"
tempname VBf
matrix `VBf' = r(V_b)
local m = mreldif(`VBf', `V4')
_cef_chk "r(V_b) = stcox efron vce(cluster id) e(V) (no ties)" `m' 1e-8

local t1 = clock("`c(current_time)'", "hms")
di _n as res "stpaf_cert_efron.do: all tests passed" as txt " (" (`t1' - `t0')/1000 " s)"
