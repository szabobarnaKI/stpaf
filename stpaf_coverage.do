* stpaf_coverage.do  --  confidence-interval coverage simulation for stpaf 0.1.1
*
* Data-generating process (same as stpaf_cert.do):
*   recurrent events, Poisson process given a gamma frailty (mean 1, var 0.5),
*   two strata (baseline rate 0.6 / 1.0), age ~ N(0,1) with log-rate 0.3,
*   time-varying exposure z(t) = 1(t > sw), sw ~ Exp(mean 2), log-rate 0.7,
*   administrative censoring C ~ U(2,5) independent of everything.
*   Because censoring does not depend on the frailty, the marginal rate is exactly
*   proportional (the LWYY model is correctly specified) and the true PAF has a
*   closed form that depends only on sw, C and tau:
*       PAF(tau) = 1 - E[min(C,tau)] / ( E[l0] + exp(0.7) E[l1] )
*   with l1 = max(min(C,tau) - sw, 0), l0 = min(C,tau) - l1.
*
* For each sample size: reps replications, model stcox z age, strata(stratum)
* vce(cluster id) breslow; stpaf, eliminate(z) horizon(2 100)  (100 = end of follow-up).
* Reported: bias, empirical SD vs mean SE, coverage of the log(1-PAF) CI (default)
* and of the untransformed CI (citrans(none)). With 1000 replications the Monte
* Carlo standard error of a 95% coverage is about 0.007.

version 15.0
clear all
set more off
set type double

local reps  1000
local sizes 50 100 400 1000

* ------------------------------------------------------------- true values
local e7 = exp(0.7)
local El1 = 3.5 - 2 + 2*(2/3)*(exp(-1) - exp(-2.5))
local true100 = 1 - 3.5/((3.5 - `El1') + `e7'*`El1')
local El1 = 2*exp(-1)
local true2 = 1 - 2/((2 - `El1') + `e7'*`El1')
di as txt "true PAF(2) = " %8.6f `true2' "   true PAF(end of follow-up) = " %8.6f `true100'

* ------------------------------------------------------------- one data set
capture program drop _cov_sim
program define _cov_sim
    syntax, n(integer)
    clear
    qui set obs `n'
    qui gen long id = _n
    qui gen double age = rnormal()
    qui gen byte stratum = 1 + (runiform() > .5)
    qui gen double C = 2 + 3*runiform()
    qui gen double sw = rexponential(2)
    qui gen double fr = rgamma(2, .5)
    qui gen double len0 = min(sw, C)
    qui gen double len1 = max(C - sw, 0)
    qui gen double r0 = (0.6 + 0.4*(stratum == 2))*fr*exp(0.3*age)
    qui gen long n0 = rpoisson(r0*len0)
    qui gen long n1 = rpoisson(r0*exp(0.7)*len1)
    qui replace n0 = 0 if n0 >= .
    qui replace n1 = 0 if n1 >= .
    qui expand n0 + n1 + 1
    qui bysort id: gen long k = _n
    qui gen double tev = .
    qui replace tev = runiform()*len0 if k <= n0
    qui replace tev = sw + runiform()*len1 if k > n0 & k <= n0 + n1
    qui replace tev = C if k == n0 + n1 + 1
    qui gen byte ev = k <= n0 + n1
    sort id tev
    qui by id: gen double tstart = cond(_n == 1, 0, tev[_n-1])
    qui drop if tev <= tstart
    qui stset tev, id(id) failure(ev) time0(tstart) exit(time .)
    qui stsplit blk, after(time = sw) at(0)
    qui gen byte z = (_t0 >= sw)
end

* ------------------------------------------------------------- simulation
tempname pf
tempfile res
postfile `pf' int n int rep double(p2 s2 l2 u2 p s l u) byte ok double sec using `res', replace
set seed 20261102
timer clear
foreach n of local sizes {
    di as txt _n "n = `n': " _c
    forvalues r = 1/`reps' {
        timer clear 1
        timer on 1
        capture {
            _cov_sim, n(`n')
            stcox z age, strata(stratum) vce(cluster id) breslow
            stpaf, eliminate(z) horizon(2 100)
        }
        local rc = _rc
        timer off 1
        if !`rc' matrix T = r(table)     // grab before timer list overwrites r()
        qui timer list 1
        local sec = r(t1)
        if `rc' {
            post `pf' (`n') (`r') (.) (.) (.) (.) (.) (.) (.) (.) (0) (`sec')
        }
        else {
            post `pf' (`n') (`r') (T[1,3]) (T[1,4]) (T[1,5]) (T[1,6]) ///
                (T[2,3]) (T[2,4]) (T[2,5]) (T[2,6]) (1) (`sec')
        }
        if mod(`r', 100) == 0 di as txt "`r' " _c
    }
}
postclose `pf'

* ------------------------------------------------------------- summary
use `res', clear
local z = invnormal(.975)
foreach t in 2 100 {
    local sfx = cond(`t' == 2, "2", "")
    gen byte cl`t' = inrange(`true`t'', l`sfx', u`sfx') if ok
    gen byte cn`t' = inrange(`true`t'', p`sfx' - `z'*s`sfx', p`sfx' + `z'*s`sfx') if ok
}
di _n as txt "{hline 92}"
di as txt "Horizon      n   fails    bias    emp.SD   mean SE  SE/SD   cover(log)  cover(none)  sec/rep"
di as txt "{hline 92}"
foreach t in 2 100 {
    local sfx = cond(`t' == 2, "2", "")
    local lab = cond(`t' == 2, "tau=2", "end")
    foreach n of local sizes {
        qui count if n == `n' & !ok
        local nf = r(N)
        qui summarize p`sfx' if n == `n' & ok
        local bias = r(mean) - `true`t''
        local sd = r(sd)
        qui summarize s`sfx' if n == `n' & ok, meanonly
        local mse = r(mean)
        qui summarize cl`t' if n == `n' & ok, meanonly
        local cl = r(mean)
        qui summarize cn`t' if n == `n' & ok, meanonly
        local cn = r(mean)
        qui summarize sec if n == `n', meanonly
        local sec = r(mean)
        di as txt %-8s "`lab'" as res %6.0f `n' %7.0f `nf' %9.4f `bias' %9.4f `sd' %9.4f `mse' ///
            %7.3f `mse'/`sd' %11.3f `cl' %12.3f `cn' %10.3f `sec'
    }
}
di as txt "{hline 92}"
di as txt "true PAF(2) = " %6.4f `true2' ", true PAF(end) = " %6.4f `true100' ///
    "; Monte Carlo SE of coverage ~ " %5.3f sqrt(.95*.05/`reps')
