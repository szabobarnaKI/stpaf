*! stpaf 0.5.1  03nov2026
*! Population attributable fraction after stcox (Cox / LWYY) with time-varying covariates.
*! Influence-function (sandwich) confidence intervals, or an exact cluster jackknife (vce(jackknife)).
*! 0.3.0: cluster() option; error for ghoshlin data without clusters
*! 0.3.1: vce(jackknife) failed for unstratified models (_stpaf_recode returned a row vector)
*! 0.4.0: fast jackknife: leave-one-cluster-out Breslow fits by Newton-Raphson in Mata (warm start
*!        from e(b), switch off with nowarm); options refit (stcox refits, as before) and nowarm
*! 0.5.0: vce(jackknife) is the default with fewer than 100 clusters (influence function if it
*!        fails); clusters whose Newton fit fails are refitted by stcox instead of stopping;
*!        weighted models with efron/exact ties use stcox refits; refit/nowarm imply vce(jackknife)
*! 0.5.1: a failing default jackknife no longer prints the Mata error trace (only the note)

program define stpaf, rclass
    version 15.0
    syntax [if] [in] [, ELIMinate(varlist numeric) SET(string) CAP(string) ///
        HORizon(numlist >0 sort) Level(cilevel) CITrans(string) EACH VCE(string) NODOTS ///
        CLuster(varname) REFIT noWARM ]

    // ------------------------------------------------------------------ checks
    st_is 2 analysis
    // stcox stores e(cmd) = "cox" and e(cmd2) = "stcox"
    if "`e(cmd2)'" != "stcox" {
        di as err "stpaf requires {bf:stcox} as the last estimation command"
        exit 301
    }
    if "`e(offset)'" != "" | "`e(shared)'" != "" {
        di as err "stcox models with offset() or shared() are not supported by stpaf"
        exit 198
    }
    local ceq : coleq e(b)
    if strpos(" `ceq' ", "tvc") {
        di as err "models fitted with tvc() are not supported; use stsplit and explicit interactions"
        exit 198
    }
    if "`citrans'" == "" local citrans "log"
    if !inlist("`citrans'", "log", "none") {
        di as err "citrans() must be {bf:log} or {bf:none}"
        exit 198
    }
    local logci = ("`citrans'" == "log")
    // variance: influence function (default) or exact cluster jackknife
    // no vce(): decided once the clusters are known (jackknife with fewer than 100 clusters)
    local vce = strtrim("`vce'")
    local vcedef 0
    if "`vce'" == "" {
        local vcedef 1
        local vce "if"
    }
    else if inlist("`vce'", "if", "influence", "robust") local vce "if"
    else if inlist("`vce'", "jackknife", "jack", "jackk", "jackkn", "jackkni", "jackknif") local vce "jackknife"
    else {
        di as err "vce() must be {bf:if} (influence function) or {bf:jackknife}"
        exit 198
    }
    // refit and nowarm are jackknife options: without vce() they request the jackknife
    if "`refit'" != "" | "`warm'" != "" {
        if `vcedef' {
            local vce "jackknife"
            local vcedef 0
        }
        else if "`vce'" != "jackknife" {
            di as err "options {bf:refit} and {bf:nowarm} require {bf:vce(jackknife)}"
            exit 198
        }
    }
    if "`refit'" != "" & "`warm'" != "" {
        di as err "options {bf:refit} and {bf:nowarm} may not be combined ({bf:stcox} refits always start from zero)"
        exit 198
    }
    if "`horizon'" == "" local horizon "."
    local nt : word count `horizon'
    local tauvec = subinstr("`horizon'", " ", ",", .)

    // ------------------------------------------------- counterfactual scenario
    local nspec 0
    foreach v of local eliminate {
        local ++nspec
        local sv`nspec' "`v'"
        local sc`nspec' "qui replace `v' = 0"
        local sl`nspec' "`v'=0"
    }
    foreach tok of local set {
        if !regexm("`tok'", "^([A-Za-z_][A-Za-z0-9_]*)=(.+)$") {
            di as err "set(): invalid specification {bf:`tok'}; use varname=# without spaces"
            exit 198
        }
        local v = regexs(1)
        local val = regexs(2)
        confirm numeric variable `v'
        confirm number `val'
        local ++nspec
        local sv`nspec' "`v'"
        local sc`nspec' "qui replace `v' = `val'"
        local sl`nspec' "`v'=`val'"
    }
    foreach tok of local cap {
        if !regexm("`tok'", "^([A-Za-z_][A-Za-z0-9_]*)=(.+)$") {
            di as err "cap(): invalid specification {bf:`tok'}; use varname=# without spaces"
            exit 198
        }
        local v = regexs(1)
        local val = regexs(2)
        confirm numeric variable `v'
        confirm number `val'
        local ++nspec
        local sv`nspec' "`v'"
        local sc`nspec' "qui replace `v' = min(`v', `val')"
        local sl`nspec' "`v'<=`val'"
    }
    if `nspec' == 0 {
        di as err "specify at least one of eliminate(), set() or cap()"
        exit 198
    }
    local nsc 1
    local jl ""
    forvalues i = 1/`nspec' {
        local scl1 `scl1' `i'
        local jl `jl' `sl`i''
    }
    local lab1 = cond(`nspec' == 1, "`sl1'", "joint")
    if "`each'" != "" & `nspec' > 1 {
        forvalues i = 1/`nspec' {
            local ++nsc
            local scl`nsc' "`i'"
            local lab`nsc' "`sl`i''"
        }
    }

    // ------------------------------------------------------- samples, ids, strata
    tempvar esamp clus strat xbchk wvar
    qui gen byte `esamp' = e(sample)
    marksample touse, novarlist
    qui replace `touse' = 0 if `esamp' != 1
    qui count if `touse'
    if r(N) == 0 error 2000

    // clusters: option cluster() > e(clustvar) > stset id() > single records.
    // The variance is built here from the Breslow information and the cluster
    // scores, so the vce() used in stcox matters only through e(clustvar).
    local clv "`cluster'"
    local clsrc "cluster()"
    if "`clv'" == "" {
        local clv "`e(clustvar)'"
        local clsrc "e(clustvar)"
    }
    if "`clv'" == "" {
        local clv : char _dta[st_id]
        local clsrc "stset id()"
    }
    local stid : char _dta[st_id]
    local stwv : char _dta[st_wv]
    if "`clv'" == "" {
        local clsrc "records"
        // ghoshlin output is stset [pw=weight_c] without id(): the records of a
        // subject would be treated as independent
        if "`stwv'" == "weight_c" {
            di as err "the data appear to be prepared by {bf:ghoshlin} (stset [pw=weight_c] without id()),"
            di as err "but no clusters are defined: specify {bf:cluster(}{it:id}{bf:)} (or {bf:cluster(repid)}"
            di as err "after expanded frequency weights), or fit {bf:stcox} with {bf:vce(cluster} {it:id}{bf:)}"
            exit 459
        }
        qui gen long `clus' = sum(`esamp') if `esamp'
    }
    else {
        capture assert !missing(`clv') if `esamp'
        if _rc {
            di as err "cluster variable {bf:`clv'} has missing values in the estimation sample"
            exit 459
        }
        qui egen long `clus' = group(`clv') if `esamp'
    }
    // cluster() should not split subjects (stset id) or the clusters of stcox
    if "`cluster'" != "" {
        local chkv "`e(clustvar)' `stid'"
        local chkv : list uniq chkv
        local chkv : list chkv - cluster
        foreach v of local chkv {
            tempvar cmin cmax
            qui egen long `cmin' = min(`clus') if `esamp', by(`v')
            qui egen long `cmax' = max(`clus') if `esamp', by(`v')
            qui count if `esamp' & `cmin' != `cmax'
            if r(N) > 0 {
                di as txt "note: {bf:`v'} is not nested within cluster({bf:`cluster'});" ///
                    " the variance is clustered on {bf:`cluster'} only"
            }
            drop `cmin' `cmax'
        }
    }
    local strv "`e(strata)'"
    if "`strv'" == "" {
        qui gen long `strat' = 1 if `esamp'
    }
    else {
        qui egen long `strat' = group(`strv') if `esamp'
    }

    // weights set with stset (fweight, pweight, iweight); they may vary across records
    local wtype "`e(wtype)'"
    if "`wtype'" == "" {
        qui gen double `wvar' = 1 if `esamp'
    }
    else {
        local wexp "`e(wexp)'"
        if `"`wexp'"' == "" {
            local wv : char _dta[st_wv]
            local wexp "= `wv'"
        }
        qui gen double `wvar' `wexp' if `esamp'
        qui count if `wvar' < 0 & `esamp'
        if r(N) > 0 {
            di as err "negative weights are not allowed"
            exit 402
        }
    }
    // frequency weights with one cluster per record: each record stands for w independent units
    local fwrep = ("`wtype'" == "fweight" & "`clv'" == "")
    // default variance: cluster jackknife with fewer than 100 clusters (the influence-function
    // CI is too narrow there), otherwise the influence function
    local jkdefault 0
    if `vcedef' {
        qui summarize `clus' if `esamp', meanonly
        if r(max) < 100 & !`fwrep' {
            local vce "jackknife"
            local jkdefault 1
        }
    }
    if `fwrep' & "`vce'" == "jackknife" {
        di as err "vce(jackknife) is not available with fweights and no clustering; use vce(if)"
        exit 198
    }
    if "`wtype'" != "" & "`clv'" == "" & "`wtype'" != "fweight" {
        di as txt "note: weighted model without clusters (no cluster(), vce(cluster) or stset id()):" ///
            " every record is treated as its own cluster."
        di as txt "      If several records belong to one subject, specify cluster({it:id})."
    }

    // ------------------------------------------------------------ design matrix
    tempname bmat omit bk
    matrix `bmat' = e(b)
    _ms_omit_info `bmat'
    matrix `omit' = r(omit)
    local cn : colnames `bmat'
    local pall = colsof(`bmat')
    local fvl ""
    local keepj ""
    local j 0
    foreach t of local cn {
        local ++j
        if `omit'[1, `j'] == 0 {
            local fvl `fvl' `t'
            local keepj `keepj' `j'
        }
    }
    local p : word count `keepj'
    if `p' == 0 {
        di as err "no non-omitted coefficients in e(b)"
        exit 322
    }
    matrix `bk' = J(1, `p', 0)
    local i 0
    foreach j of local keepj {
        local ++i
        matrix `bk'[1, `i'] = `bmat'[1, `j']
    }
    fvrevar `fvl'
    local xvars "`r(varlist)'"
    qui predict double `xbchk' if `esamp', xb

    // ------------------------------- counterfactual design matrices (Mata globals)
    // Each scenario's design matrix is rebuilt in Stata (so factor variables and
    // interactions are recomputed) and stored as a plain Mata matrix
    // __stpaf_Xc1, __stpaf_Xc2, ...; all computations then run in one Mata call.
    tempname TAB SC ST CHK
    forvalues j = 1/`nsc' {
        capture mata: mata drop __stpaf_Xc`j'
    }
    capture noisily {
        forvalues j = 1/`nsc' {
            preserve
            foreach i of local scl`j' {
                capture recast double `sv`i''
                `sc`i''
            }
            fvrevar `fvl'
            local xcv "`r(varlist)'"
            mata: __stpaf_Xc`j' = st_data(., tokens("`xcv'"), "`esamp'")
            restore
        }
        mata: _stpaf_run(st_data(., tokens("`xvars'"), "`esamp'"), ///
            st_data(., "_t0", "`esamp'"), st_data(., "_t", "`esamp'"), ///
            st_data(., "_d", "`esamp'"), st_data(., "`touse'", "`esamp'"), ///
            st_data(., "`wvar'", "`esamp'"), ///
            st_data(., "`clus'", "`esamp'"), st_data(., "`strat'", "`esamp'"), ///
            st_matrix("`bk'"), `nsc', (`tauvec'), `level', `logci', `fwrep', ///
            "`TAB'", "`SC'", "`ST'", "`CHK'", "`xbchk'", "`esamp'")
    }
    local rc = _rc
    if !`rc' & scalar(`CHK') > 1e-6 {
        di as err "stpaf could not reproduce the linear predictor of stcox (max abs diff = " ///
            %9.3g scalar(`CHK') ")"
        di as err "the model probably contains terms that cannot be rebuilt from e(b)"
        local rc 498
    }

    // ----------------------------------------------------- exact cluster jackknife
    // Each cluster is left out in turn and the model is re-estimated on the reduced
    // sample (same terms, strata and stset weights); the PAFs of the reduced samples
    // give the jackknife variance of log(1-PAF) and of PAF, which replaces the
    // influence-function one.
    // Default (Breslow ties, or an unweighted model without tied failure times):
    // Newton-Raphson for the Breslow partial likelihood in Mata, the left-out cluster given
    // zero weight, started from e(b) (nowarm: from 0) and iterated until mreldif(b) <= 1e-12
    // and U'I^-1U <= 1e-12, i.e. to convergence independent of the starting values.
    // Clusters whose Newton fit fails (singular information, no convergence) are refitted
    // by stcox. Option refit, efron/exact ties with tied failure times, and weighted
    // models with efron/exact ties: stcox refits for all clusters.
    // When the jackknife was chosen by default (fewer than 100 clusters) and fails, the
    // influence-function results are kept.
    local jkmeth ""
    local jkit = .
    local jknref 0
    local jkfail 0
    if !`rc' & "`vce'" == "jackknife" {
        local G0 = `ST'[1, 1]
        local jkmeth "refit"
        local full 1
        tempname JKI JKNF JKF
        capture mata: mata drop __stpaf_PF
        // a jackknife chosen by default falls back to the influence function on error, so its
        // Mata error messages are not shown (noisily would also override quietly)
        local capn = cond(`jkdefault', "capture", "capture noisily")
        if "`refit'" == "" {
            local brestie = inlist("`e(ties)'", "breslow", "")
            local dowarm = ("`warm'" == "")
            if !`brestie' & "`wtype'" != "" {
                // defensive: stcox 19 allows no weights with efron (cert 12d); if a weighted
                // efron/exact fit occurs, fweights > 1 imply tied failures that a count of
                // failure records cannot see, so refit with the model's own ties
                di as txt "note: weighted model with {bf:`e(ties)'} ties; the jackknife refits {bf:stcox}"
                if "`warm'" != "" di as txt "      (option nowarm has no effect)"
            }
            else {
                `capn' mata: _stpaf_jackfast(st_data(., tokens("`xvars'"), "`esamp'"), ///
                    st_data(., "_t0", "`esamp'"), st_data(., "_t", "`esamp'"), ///
                    st_data(., "_d", "`esamp'"), st_data(., "`touse'", "`esamp'"), ///
                    st_data(., "`wvar'", "`esamp'"), ///
                    st_data(., "`clus'", "`esamp'"), st_data(., "`strat'", "`esamp'"), ///
                    st_matrix("`bk'"), `nsc', (`tauvec'), `level', `logci', "`TAB'", ///
                    `dowarm', `brestie', "`JKI'", "`JKNF'", "`JKF'")
                local rc = _rc
                if !`rc' {
                    if scalar(`JKI') >= 0 {
                        local jkmeth = cond(`dowarm', "newton", "newton, nowarm")
                        local jkit = scalar(`JKI')
                        local jknref = scalar(`JKNF')
                        local full 0
                        if `jknref' > 0 {
                            di as txt "note: the Newton fit failed for `jknref' of `G0' left-out clusters" ///
                                " (singular information or no convergence);"
                            di as txt "      {bf:stcox} is refitted for " cond(`jknref' == 1, "that cluster", "these clusters")
                        }
                    }
                    else {
                        di as txt "note: tied failure times with {bf:`e(ties)'} ties; the jackknife refits {bf:stcox}"
                        if "`warm'" != "" di as txt "      (option nowarm has no effect)"
                    }
                }
            }
        }
        // stcox refits: all clusters (full) or those whose Newton fit failed (JKF)
        local nref = cond(`full', `G0', `jknref')
        if !`rc' & `nref' > 0 {
            local ties "`e(ties)'"
            if !inlist("`ties'", "breslow", "efron", "exactm", "exactp") local ties "breslow"
            local stopt = cond("`strv'" == "", "", "strata(`strv')")
            tempname BJ H
            matrix `BJ' = J(`nref', `p', .)
            _estimates hold `H', restore
            if "`nodots'" == "" {
                _dots 0, title(Jackknife replications) reps(`nref')
            }
            forvalues i = 1/`nref' {
                if `full' local g `i'
                else local g = `JKF'[`i', 1]
                capture quietly stcox `fvl' if `esamp' & `clus' != `g', `stopt' `ties'
                local rcg = _rc
                if !`rcg' {
                    local j 0
                    foreach t of local fvl {
                        local ++j
                        matrix `BJ'[`i', `j'] = _b[`t']
                    }
                }
                if "`nodots'" == "" _dots `i' `rcg'
                if `rcg' {
                    di
                    di as err "stcox failed when cluster `g' of `G0' was left out (rc = `rcg')"
                    local rc `rcg'
                    continue, break
                }
            }
            if "`nodots'" == "" di
            _estimates unhold `H'
            if !`rc' {
                if `full' {
                    `capn' mata: _stpaf_jack(st_data(., tokens("`xvars'"), "`esamp'"), ///
                        st_data(., "_t0", "`esamp'"), st_data(., "_t", "`esamp'"), ///
                        st_data(., "_d", "`esamp'"), st_data(., "`touse'", "`esamp'"), ///
                        st_data(., "`wvar'", "`esamp'"), ///
                        st_data(., "`clus'", "`esamp'"), st_data(., "`strat'", "`esamp'"), ///
                        st_matrix("`BJ'"), `nsc', (`tauvec'), `level', `logci', "`TAB'")
                }
                else {
                    `capn' mata: _stpaf_jackfill(st_data(., tokens("`xvars'"), "`esamp'"), ///
                        st_data(., "_t0", "`esamp'"), st_data(., "_t", "`esamp'"), ///
                        st_data(., "_d", "`esamp'"), st_data(., "`touse'", "`esamp'"), ///
                        st_data(., "`wvar'", "`esamp'"), ///
                        st_data(., "`clus'", "`esamp'"), st_data(., "`strat'", "`esamp'"), ///
                        st_matrix("`BJ'"), st_matrix("`JKF'"), `nsc', (`tauvec'), `level', `logci', "`TAB'")
                }
                local rc = _rc
            }
        }
        capture mata: mata drop __stpaf_PF
        // a jackknife chosen by default must not stop stpaf: keep the influence function
        if `rc' & `rc' != 1 & `jkdefault' {
            di as txt "note: the cluster jackknife (default with fewer than 100 clusters) failed (rc = `rc');"
            di as txt "      standard errors from the influence function are reported instead"
            local rc 0
            local vce "if"
            local jkfail 1
            local jkmeth ""
            local jknref 0
        }
    }
    forvalues j = 1/`nsc' {
        capture mata: mata drop __stpaf_Xc`j'
    }
    if `rc' exit `rc'

    // ---------------------------------------------------------------- display
    local G    = `ST'[1, 1]
    local nev  = `ST'[1, 3]
    local nK   = `ST'[1, 4]
    local nS   = `ST'[1, 5]
    local maxsc = 0
    forvalues j = 1/`p' {
        if `SC'[1, `j'] > `maxsc' local maxsc = `SC'[1, `j']
    }

    di
    di as txt "Population attributable fraction after {bf:stcox}" ///
        _col(55) "Clusters" _col(68) "= " as res %9.0fc `G'
    di as txt "Estimand: 1 - E[cf. events] / E[observed events]" ///
        _col(55) "Failures" _col(68) "= " as res %9.0fc `nev'
    di as txt "Counterfactual: " as res abbrev("`jl'", 34) ///
        as txt _col(55) "Event times" _col(68) "= " as res %9.0fc `nK'
    if `nS' > 1 {
        di as txt "Baseline hazards stratified by: " as res abbrev("`strv'", 22) ///
            as txt _col(55) "Strata" _col(68) "= " as res %9.0fc `nS'
    }
    if "`wtype'" != "" {
        di as txt "Weights: " as res "[`wtype' `wexp']" as txt "   (failures are weighted counts)"
    }
    if "`vce'" == "jackknife" {
        if "`jkmeth'" == "refit" local jkdesc "`G' stcox refits"
        else {
            local jkdesc = "`=`G'-`jknref'' Newton fits, " + cond("`warm'" == "", "warm start", "cold start") + ///
                ", mean " + string(`jkit', "%4.2f") + " iterations"
            if `jknref' > 0 local jkdesc = "`jkdesc', `jknref' stcox refit" + cond(`jknref' > 1, "s", "")
        }
        di as txt "Std. err. from the cluster jackknife (" as res "`jkdesc'" as txt "), clusters: " ///
            as res cond("`clv'" == "", "record", "`clv'") as txt " (`clsrc')"
        if `jkdefault' {
            di as txt "(the jackknife is the default with fewer than 100 clusters; {bf:vce(if)} gives the influence function)"
        }
    }
    else {
        di as txt "Std. err. from the influence function, clustered on " ///
            as res cond("`clv'" == "", "record", "`clv'") as txt " (`clsrc')"
    }
    di
    di as txt "{hline 13}{c TT}{hline 70}"
    di as txt "{ralign 12:Scenario} {c |}" _col(17) "Horizon" _col(28) "PAF" _col(37) "Std.Err." ///
        _col(48) "[`level'% Conf. Interval]" _col(72) "E[obs]" _col(82) "E[cf]"
    di as txt "{hline 13}{c +}{hline 70}"
    local nr = rowsof(`TAB')
    forvalues i = 1/`nr' {
        local sj = `TAB'[`i', 1]
        local h  = `TAB'[`i', 2]
        di as txt "{ralign 12:" abbrev("`lab`sj''", 12) "} {c |}" as res _col(17) %7.4g `h' ///
            _col(26) %8.4f `TAB'[`i', 3] _col(36) %8.4f `TAB'[`i', 4] ///
            _col(48) "[" %7.4f `TAB'[`i', 5] "," %7.4f `TAB'[`i', 6] "]" ///
            _col(70) %8.1f `TAB'[`i', 7] _col(80) %8.1f `TAB'[`i', 8]
    }
    di as txt "{hline 13}{c BT}{hline 70}"
    if "`citrans'" == "log" {
        di as txt "CI is computed on the log(1-PAF) scale and back-transformed."
    }
    if "`vce'" == "if" & `G' < 100 & !`jkfail' & !`fwrep' {
        di as txt "Note: with fewer than 100 clusters the influence-function CI tends to be too narrow;"
        di as txt "      {bf:vce(jackknife)}, the default in this case, is recommended."
    }
    if `maxsc' > 1e-4 & "`vce'" == "if" {
        di as txt "Note: standardized partial-likelihood score at e(b) = " %9.2g `maxsc' ///
            "; the model was probably fitted with efron or exactp ties."
        di as txt "      Breslow-type expressions are used for the variance; refit with the" ///
            " breslow option for exact agreement."
    }

    // --------------------------------------------------------------- returns
    local last = `nt'
    return scalar paf    = `TAB'[`last', 3]
    return scalar se     = `TAB'[`last', 4]
    return scalar lb     = `TAB'[`last', 5]
    return scalar ub     = `TAB'[`last', 6]
    return scalar E      = `TAB'[`last', 7]
    return scalar E_cf   = `TAB'[`last', 8]
    return scalar se_logR = `TAB'[`last', 9]
    return scalar N_clust = `G'
    return scalar N_fail  = `nev'
    return scalar level   = `level'
    matrix colnames `TAB' = scenario horizon paf se lb ub E E_cf se_logR events
    return matrix table = `TAB'
    return local cfactual "`jl'"
    return local citrans "`citrans'"
    return local vce "`vce'"
    return local clustvar "`clv'"
    if "`vce'" == "jackknife" {
        return local jk_method "`jkmeth'"
        if "`jkmeth'" != "refit" {
            return scalar jk_iter = `jkit'
            return scalar jk_nrefit = `jknref'
        }
    }
    return local cmd "stpaf"
end

// =============================================================================
// Mata
// =============================================================================
mata:

struct stpaf_dat {
    real scalar    n
    real scalar    p
    real scalar    G
    real scalar    S
    real scalar    Ktot
    real scalar    cen
    real matrix    Eb
    real matrix    Gc
    real matrix    Dc
    real matrix    B
    real matrix    IFb
    real matrix    info
    real colvector xb
    real colvector a
    real colvector w
    real colvector t0
    real colvector t1
    real colvector d
    real colvector g
    real colvector clus
    real colvector strat
    real colvector ord
    real colvector tk
    real colvector dk
    real colvector Kst
    real colvector gbase
    real colvector S0k
    real colvector hk
    real colvector Lc
    real colvector Hc
    real colvector DHc
    real colvector DKc
    real colvector rm0
    real colvector rm1
    real rowvector b
    real rowvector scorechk
}

// number of elements of tk[off+1..off+K] (sorted ascending) that are <= t
real scalar _stpaf_pos(real colvector tk, real scalar off, real scalar K, real scalar t)
{
    real scalar lo, hi, mid
    lo = 0
    hi = K
    while (lo < hi) {
        mid = floor((lo + hi + 1)/2)
        if (tk[off + mid] <= t) lo = mid
        else hi = mid - 1
    }
    return(lo)
}

struct stpaf_dat scalar _stpaf_prep(real matrix X, real colvector t0, real colvector t1,
    real colvector d, real colvector g, real colvector w, real colvector clus, real colvector strat,
    real rowvector b)
{
    struct stpaf_dat scalar D
    real scalar n, p, G, S, Ktot, nE, i, s, k, r, ns, Ks, gb, pb, ig, i1, i0, m0, m1
    real scalar S0, T0, tkk, ar, lam
    real colvector eta, a, ev, es, et, ew, oE, kid, newk, ks, ord, rs, o1, o0, tt1, tt0
    real colvector tk, dk, Kst, gbase, S0k, T0k, hk, Lc, Hc, DHc, DKc, rm0, rm1
    real rowvector xr, S1, T1, ur
    real matrix S2, I, info, Eb, Ak, Gc, Dc, Ucl

    n = rows(X)
    p = cols(X)
    G = max(clus)
    S = max(strat)
    if (min(clus) < 1 | min(strat) < 1) _error("invalid cluster or strata codes")
    if (G < 2) _error("at least 2 clusters are required")

    eta = X*b'
    D.xb = eta
    D.cen = mean(eta)
    if (min(w) < 0) _error("negative weights")
    a = w :* exp(eta :- D.cen)

    // ---- event-time grid, by stratum
    ev = selectindex((d :== 1) :& (w :> 0))
    nE = rows(ev)
    if (nE == 0) _error("no failure events in the estimation sample")
    es = strat[ev]
    et = t1[ev]
    ew = w[ev]
    oE = order((es, et), (1, 2))
    es = es[oE]
    et = et[oE]
    ew = ew[oE]
    newk = J(nE, 1, 1)
    if (nE > 1) {
        newk[|2 \ nE|] = (es[|2 \ nE|] :!= es[|1 \ nE-1|]) :| (et[|2 \ nE|] :!= et[|1 \ nE-1|])
    }
    kid = runningsum(newk)
    Ktot = kid[nE]
    tk = J(Ktot, 1, .)
    dk = J(Ktot, 1, 0)
    ks = J(Ktot, 1, .)
    for (i = 1; i <= nE; i++) {
        tk[kid[i]] = et[i]
        ks[kid[i]] = es[i]
        dk[kid[i]] = dk[kid[i]] + ew[i]
    }
    Kst = J(S, 1, 0)
    for (k = 1; k <= Ktot; k++) Kst[ks[k]] = Kst[ks[k]] + 1
    gbase = J(S, 1, 0)
    for (s = 2; s <= S; s++) gbase[s] = gbase[s-1] + Kst[s-1]

    ord = order(strat, 1)
    info = panelsetup(strat[ord], 1)
    if (rows(info) != S) _error("strata coding error")

    // ---- descending sweep over event times: risk-set sums at every event time
    S0k = J(Ktot, 1, 0)
    T0k = J(Ktot, 1, 0)
    hk = J(Ktot, 1, 0)
    Eb = J(Ktot, p, 0)
    Ak = J(Ktot, p, 0)
    I = J(p, p, 0)
    rm0 = J(n, 1, 0)
    rm1 = J(n, 1, 0)
    for (s = 1; s <= S; s++) {
        Ks = Kst[s]
        if (Ks == 0) continue
        gb = gbase[s]
        rs = ord[|info[s,1] \ info[s,2]|]
        ns = rows(rs)
        for (i = 1; i <= ns; i++) {
            r = rs[i]
            rm0[r] = _stpaf_pos(tk, gb, Ks, t0[r])
            rm1[r] = _stpaf_pos(tk, gb, Ks, t1[r])
        }
        tt1 = t1[rs]
        tt0 = t0[rs]
        o1 = order(tt1, -1)
        o0 = order(tt0, -1)
        S0 = 0
        T0 = 0
        S1 = J(1, p, 0)
        T1 = J(1, p, 0)
        S2 = J(p, p, 0)
        i1 = 1
        i0 = 1
        for (k = Ks; k >= 1; k--) {
            tkk = tk[gb + k]
            while (i1 <= ns) {
                r = rs[o1[i1]]
                if (t1[r] < tkk) break
                ar = a[r]
                xr = X[r,.]
                S0 = S0 + ar
                S1 = S1 + ar*xr
                S2 = S2 + ar*(xr'*xr)
                if (g[r]) {
                    T0 = T0 + ar
                    T1 = T1 + ar*xr
                }
                i1++
            }
            while (i0 <= ns) {
                r = rs[o0[i0]]
                if (t0[r] < tkk) break
                ar = a[r]
                xr = X[r,.]
                S0 = S0 - ar
                S1 = S1 - ar*xr
                S2 = S2 - ar*(xr'*xr)
                if (g[r]) {
                    T0 = T0 - ar
                    T1 = T1 - ar*xr
                }
                i0++
            }
            ig = gb + k
            S0k[ig] = S0
            T0k[ig] = T0
            Eb[ig,.] = S1/S0
            Ak[ig,.] = T1/S0
            I = I + dk[ig]*(S2/S0 - (S1'*S1)/(S0*S0))
        }
    }

    // ---- prefix (cumulative) arrays, one block of length K_s+1 per stratum;
    //      element m (0..K_s) of stratum s lives at index gbase[s] + s + m
    Lc = J(Ktot + S, 1, 0)
    Hc = J(Ktot + S, 1, 0)
    DHc = J(Ktot + S, 1, 0)
    DKc = J(Ktot + S, 1, 0)
    Gc = J(Ktot + S, p, 0)
    Dc = J(Ktot + S, p, 0)
    for (s = 1; s <= S; s++) {
        Ks = Kst[s]
        gb = gbase[s]
        pb = gb + s
        for (k = 1; k <= Ks; k++) {
            ig = gb + k
            lam = dk[ig]/S0k[ig]
            hk[ig] = T0k[ig]/S0k[ig]
            Lc[pb+k] = Lc[pb+k-1] + lam
            Hc[pb+k] = Hc[pb+k-1] + hk[ig]*lam
            DHc[pb+k] = DHc[pb+k-1] + dk[ig]*hk[ig]
            DKc[pb+k] = DKc[pb+k-1] + dk[ig]
            Gc[pb+k,.] = Gc[pb+k-1,.] + lam*Eb[ig,.]
            Dc[pb+k,.] = Dc[pb+k-1,.] + dk[ig]*(Ak[ig,.] - hk[ig]*Eb[ig,.])
        }
    }

    // ---- cluster-level partial-likelihood score (Breslow)
    Ucl = J(G, p, 0)
    for (r = 1; r <= n; r++) {
        s = strat[r]
        gb = gbase[s]
        pb = gb + s
        m0 = rm0[r]
        m1 = rm1[r]
        ur = -a[r]*(X[r,.]*(Lc[pb+m1] - Lc[pb+m0]) - (Gc[pb+m1,.] - Gc[pb+m0,.]))
        if (d[r] & w[r] > 0) ur = ur + w[r]*(X[r,.] - Eb[gb+m1,.])
        Ucl[clus[r],.] = Ucl[clus[r],.] + ur
    }

    D.n = n
    D.p = p
    D.G = G
    D.S = S
    D.Ktot = Ktot
    D.B = invsym(I)
    D.IFb = Ucl*D.B
    D.scorechk = abs(colsum(Ucl)) :/ sqrt(diagonal(I))'
    D.Eb = Eb
    D.Gc = Gc
    D.Dc = Dc
    D.info = info
    D.a = a
    D.w = w
    D.t0 = t0
    D.t1 = t1
    D.d = d
    D.g = g
    D.clus = clus
    D.strat = strat
    D.ord = ord
    D.tk = tk
    D.dk = dk
    D.Kst = Kst
    D.gbase = gbase
    D.S0k = S0k
    D.hk = hk
    D.Lc = Lc
    D.Hc = Hc
    D.DHc = DHc
    D.DKc = DKc
    D.rm0 = rm0
    D.rm1 = rm1
    D.b = b
    return(D)
}

// PAF (one row per horizon) for the counterfactual design matrix Xc
// columns: horizon, paf, se, lb, ub, E, E0, se(log R), failures up to horizon
real matrix _stpaf_scen(struct stpaf_dat scalar D, real matrix Xc, real rowvector tau,
    real scalar level, real scalar logci, real scalar fwrep)
{
    real scalar n, p, G, S, Ktot, nt, s, k, r, ns, Ks, gb, pb, ig, i1, i0, m0, m1, q0, q1
    real scalar T0, tkk, lam, ti, E, E0, evn, v, v0, R, z, paf, se, selr, lb, ub, mtau, Nw, vv
    real colvector ac, T0ck, h0k, H0c, DH0c, rs, o1, o0, mt, psi, psi0, phi, ii
    real rowvector T1, Dt, D0t
    real matrix Ack, D0c, res

    n = D.n
    p = D.p
    G = D.G
    S = D.S
    Ktot = D.Ktot
    nt = cols(tau)

    ac = D.w :* D.g :* exp(Xc*D.b' :- D.cen)

    // ---- counterfactual risk-set sums at every event time
    T0ck = J(Ktot, 1, 0)
    Ack = J(Ktot, p, 0)
    for (s = 1; s <= S; s++) {
        Ks = D.Kst[s]
        if (Ks == 0) continue
        gb = D.gbase[s]
        rs = D.ord[|D.info[s,1] \ D.info[s,2]|]
        ns = rows(rs)
        o1 = order(D.t1[rs], -1)
        o0 = order(D.t0[rs], -1)
        T0 = 0
        T1 = J(1, p, 0)
        i1 = 1
        i0 = 1
        for (k = Ks; k >= 1; k--) {
            tkk = D.tk[gb + k]
            while (i1 <= ns) {
                r = rs[o1[i1]]
                if (D.t1[r] < tkk) break
                T0 = T0 + ac[r]
                T1 = T1 + ac[r]*Xc[r,.]
                i1++
            }
            while (i0 <= ns) {
                r = rs[o0[i0]]
                if (D.t0[r] < tkk) break
                T0 = T0 - ac[r]
                T1 = T1 - ac[r]*Xc[r,.]
                i0++
            }
            ig = gb + k
            T0ck[ig] = T0
            Ack[ig,.] = T1/D.S0k[ig]
        }
    }
    h0k = T0ck :/ D.S0k
    H0c = J(Ktot + S, 1, 0)
    DH0c = J(Ktot + S, 1, 0)
    D0c = J(Ktot + S, p, 0)
    for (s = 1; s <= S; s++) {
        Ks = D.Kst[s]
        gb = D.gbase[s]
        pb = gb + s
        for (k = 1; k <= Ks; k++) {
            ig = gb + k
            lam = D.dk[ig]/D.S0k[ig]
            H0c[pb+k] = H0c[pb+k-1] + h0k[ig]*lam
            DH0c[pb+k] = DH0c[pb+k-1] + D.dk[ig]*h0k[ig]
            D0c[pb+k,.] = D0c[pb+k-1,.] + D.dk[ig]*(Ack[ig,.] - h0k[ig]*D.Eb[ig,.])
        }
    }

    // ---- point estimate and influence function for each horizon
    z = invnormal(1 - (1 - level/100)/2)
    res = J(nt, 9, .)
    mt = J(S, 1, 0)
    for (ti = 1; ti <= nt; ti++) {
        E = 0
        E0 = 0
        evn = 0
        Dt = J(1, p, 0)
        D0t = J(1, p, 0)
        for (s = 1; s <= S; s++) {
            Ks = D.Kst[s]
            gb = D.gbase[s]
            pb = gb + s
            mt[s] = 0
            if (Ks > 0) mt[s] = _stpaf_pos(D.tk, gb, Ks, tau[ti])
            E = E + D.DHc[pb + mt[s]]
            E0 = E0 + DH0c[pb + mt[s]]
            evn = evn + D.DKc[pb + mt[s]]
            Dt = Dt + D.Dc[pb + mt[s],.]
            D0t = D0t + D0c[pb + mt[s],.]
        }
        if (E <= 0 | E0 <= 0) _error("no expected events in the selected subset/horizon")
        psi = J(G, 1, 0)
        psi0 = J(G, 1, 0)
        for (r = 1; r <= n; r++) {
            s = D.strat[r]
            gb = D.gbase[s]
            pb = gb + s
            mtau = mt[s]
            m0 = D.rm0[r]
            m1 = D.rm1[r]
            q0 = m0
            if (q0 > mtau) q0 = mtau
            q1 = m1
            if (q1 > mtau) q1 = mtau
            v  = -D.a[r]*(D.Hc[pb+q1] - D.Hc[pb+q0]) + D.g[r]*D.a[r]*(D.Lc[pb+q1] - D.Lc[pb+q0])
            v0 = -D.a[r]*(H0c[pb+q1] - H0c[pb+q0]) + ac[r]*(D.Lc[pb+q1] - D.Lc[pb+q0])
            if (D.d[r] & (m1 <= mtau) & D.w[r] > 0) {
                v  = v  + D.w[r]*D.hk[gb+m1]
                v0 = v0 + D.w[r]*h0k[gb+m1]
            }
            psi[D.clus[r]]  = psi[D.clus[r]]  + v
            psi0[D.clus[r]] = psi0[D.clus[r]] + v0
        }
        psi  = psi  + D.IFb*Dt'
        psi0 = psi0 + D.IFb*D0t'
        phi = psi0/E0 - psi/E
        if (fwrep) {
            ii = selectindex(D.w :> 0)
            Nw = sum(D.w)
            vv = Nw/(Nw-1)*sum(phi[ii]:^2 :/ D.w[ii])
        }
        else vv = G/(G-1)*quadcross(phi, phi)
        selr = sqrt(vv)
        R = E0/E
        paf = 1 - R
        se = R*selr
        if (logci) {
            lb = 1 - R*exp(z*selr)
            ub = 1 - R*exp(-z*selr)
        }
        else {
            lb = paf - z*se
            ub = paf + z*se
        }
        res[ti,1] = tau[ti]
        if (tau[ti] >= .) res[ti,1] = max(D.tk)
        res[ti,2] = paf
        res[ti,3] = se
        res[ti,4] = lb
        res[ti,5] = ub
        res[ti,6] = E
        res[ti,7] = E0
        res[ti,8] = selr
        res[ti,9] = evn
    }
    return(res)
}

// driver: prepare once, evaluate every scenario, post results to Stata
// (the structure never leaves Mata function scope)
void _stpaf_run(real matrix X, real colvector t0, real colvector t1, real colvector d,
    real colvector g, real colvector w, real colvector clus, real colvector strat,
    real rowvector b, real scalar nsc, real rowvector tau, real scalar level,
    real scalar logci, real scalar fwrep, string scalar tabname, string scalar scname,
    string scalar stname, string scalar chkname, string scalar xbvar, string scalar sel)
{
    struct stpaf_dat scalar D
    pointer(real matrix) scalar pX
    real matrix res, out
    real scalar j

    D = _stpaf_prep(X, t0, t1, d, g, w, clus, strat, b)
    st_numscalar(chkname, _stpaf_chk(D, xbvar, sel))
    _stpaf_post(D, scname, stname)
    out = J(0, 10, .)
    for (j = 1; j <= nsc; j++) {
        pX = findexternal("__stpaf_Xc" + strofreal(j))
        if (pX == NULL) _error("stpaf internal error: counterfactual design matrix not found")
        if (rows(*pX) != rows(X) | cols(*pX) != cols(X)) {
            _error("stpaf internal error: counterfactual design matrix has wrong dimensions")
        }
        res = _stpaf_scen(D, *pX, tau, level, logci, fwrep)
        out = out \ (J(rows(res), 1, j), res)
    }
    st_matrix(tabname, out)
}

// recode positive integer codes to 1..k (keeps order)
real colvector _stpaf_recode(real colvector v)
{
    real colvector u, map
    u = uniqrows(v)
    map = J(max(u), 1, 0)
    map[u] = (1::rows(u))
    // list subscript map[v, 1] always returns a column; map[v] would return a
    // row when map is 1 x 1 (a single stratum or cluster code)
    return(map[v, 1])
}

// refit path: replicate PAFs (one row per cluster in cl, columns in the row order of TAB)
// from the coefficients refitted by stcox without that cluster (row i of BJ for cl[i])
real matrix _stpaf_jkrefit(real matrix X, real colvector t0, real colvector t1, real colvector d,
    real colvector g, real colvector w, real colvector clus, real colvector strat,
    real matrix BJ, real colvector cl, real scalar nsc, real rowvector tau, real scalar level,
    real scalar logci)
{
    struct stpaf_dat scalar D
    pointer(real matrix) scalar pX
    real matrix PF, res
    real colvector ii
    real scalar nt, i, gg, j

    nt = cols(tau)
    PF = J(rows(cl), nsc*nt, .)
    for (i = 1; i <= rows(cl); i++) {
        gg = cl[i]
        ii = selectindex(clus :!= gg)
        D = _stpaf_prep(X[ii,.], t0[ii], t1[ii], d[ii], g[ii], w[ii],
            _stpaf_recode(clus[ii]), _stpaf_recode(strat[ii]), BJ[i,.])
        for (j = 1; j <= nsc; j++) {
            pX = findexternal("__stpaf_Xc" + strofreal(j))
            if (pX == NULL) _error("stpaf internal error: counterfactual design matrix not found")
            res = _stpaf_scen(D, (*pX)[ii,.], tau, level, logci, 0)
            PF[i, ((j-1)*nt+1)..(j*nt)] = res[.,2]'
        }
    }
    return(PF)
}

// exact cluster jackknife, refit path: BJ holds the coefficients refitted by stcox
// without cluster g (row g); the PAFs of the reduced samples replace the IF-based SEs
void _stpaf_jack(real matrix X, real colvector t0, real colvector t1, real colvector d,
    real colvector g, real colvector w, real colvector clus, real colvector strat,
    real matrix BJ, real scalar nsc, real rowvector tau, real scalar level,
    real scalar logci, string scalar tabname)
{
    real matrix PF

    PF = _stpaf_jkrefit(X, t0, t1, d, g, w, clus, strat, BJ, (1::rows(BJ)), nsc, tau,
        level, logci)
    _stpaf_jkvar(PF, level, logci, tabname)
}

// fast jackknife with failed Newton fits: the replicate PAFs of the other clusters are in
// the external __stpaf_PF; those of the clusters cl come from the stcox refits BJ
void _stpaf_jackfill(real matrix X, real colvector t0, real colvector t1, real colvector d,
    real colvector g, real colvector w, real colvector clus, real colvector strat,
    real matrix BJ, real colvector cl, real scalar nsc, real rowvector tau, real scalar level,
    real scalar logci, string scalar tabname)
{
    pointer() scalar pPF
    real matrix PF

    pPF = findexternal("__stpaf_PF")
    if (pPF == NULL) _error("stpaf internal error: replicate PAFs of the Newton jackknife not found")
    PF = *pPF
    if (rows(BJ) != rows(cl)) _error("stpaf internal error: refitted clusters do not match")
    PF[cl, .] = _stpaf_jkrefit(X, t0, t1, d, g, w, clus, strat, BJ, cl, nsc, tau, level, logci)
    _stpaf_jkvar(PF, level, logci, tabname)
}

// jackknife variance from the replicate PAFs (G x (scenarios*horizons), in the
// row order of TAB): (G-1)/G sum (theta_-g - mean)^2 of log(1-PAF) and of PAF
void _stpaf_jkvar(real matrix PF, real scalar level, real scalar logci, string scalar tabname)
{
    real matrix T
    real colvector lr, pf
    real scalar G, c, z, R, selr, sep, paf

    G = rows(PF)
    if (hasmissing(PF)) _error("jackknife replicate PAF could not be computed")
    if (max(PF) >= 1) _error("jackknife replicate PAF >= 1; log(1-PAF) undefined")
    T = st_matrix(tabname)
    z = invnormal(1 - (1 - level/100)/2)
    for (c = 1; c <= cols(PF); c++) {
        pf = PF[., c]
        lr = ln(1 :- pf)
        selr = sqrt((G-1)/G*sum((lr :- mean(lr)):^2))
        sep  = sqrt((G-1)/G*sum((pf :- mean(pf)):^2))
        paf = T[c, 3]
        R = 1 - paf
        T[c, 4] = sep
        T[c, 9] = selr
        if (logci) {
            T[c, 5] = 1 - R*exp(z*selr)
            T[c, 6] = 1 - R*exp(-z*selr)
        }
        else {
            T[c, 5] = paf - z*sep
            T[c, 6] = paf + z*sep
        }
    }
    st_matrix(tabname, T)
}

// ----------------------------------------------------------------- fast jackknife
// Leave-one-cluster-out re-estimation without refitting in Stata: the left-out cluster
// gets weight 0, so the event grid, the records' positions on it (rm0, rm1) and all
// orderings are those of the full sample and are set up once. Risk-set sums at the
// event times are bucket sums (panelsum) of the records over the event at which they
// leave the risk set (last event at risk, m1) minus those over the event before they
// enter it (m0), cumulated backwards within strata.

struct stpaf_jk {
    real scalar    n
    real scalar    p
    real scalar    S
    real scalar    Ktot
    real scalar    cen
    real colvector Kst
    real colvector gbase
    real colvector gb
    real colvector strat
    real colvector m0
    real colvector m1
    real colvector cp0
    real colvector cp1
    real colvector o1
    real colvector o0
    real colvector of
    real colvector k1
    real colvector k0
    real colvector kf
    real matrix    in1
    real matrix    in0
    real matrix    inf
    real matrix    X1
    real matrix    X0
    real matrix    Xf
}

// sort the records idx by their event index base + m; panel info and the event index
// of each panel; Xs = X of the sorted records
void _stpaf_jkpanel(real colvector idx, real colvector base, real colvector m, real matrix X,
    real colvector o, real matrix info, real colvector k, real matrix Xs)
{
    real colvector key, ord, ks

    if (rows(idx) == 0) {
        o = J(0, 1, .)
        info = J(0, 2, .)
        k = J(0, 1, .)
        Xs = J(0, cols(X), .)
        return
    }
    key = base[idx, 1] + m[idx, 1]
    ord = order(key, 1)
    o = idx[ord, 1]
    ks = key[ord, 1]
    info = panelsetup(ks, 1)
    k = ks[info[., 1], 1]
    Xs = X[o, .]
}

struct stpaf_jk scalar _stpaf_jkset(real matrix X, struct stpaf_dat scalar D)
{
    struct stpaf_jk scalar JK
    real colvector act, idx, o, k
    real matrix info, Xs

    JK.n = D.n
    JK.p = D.p
    JK.S = D.S
    JK.Ktot = D.Ktot
    JK.cen = D.cen
    JK.Kst = D.Kst
    JK.gbase = D.gbase
    JK.strat = D.strat
    JK.m0 = D.rm0
    JK.m1 = D.rm1
    // list subscripts [., 1] keep column vectors even when there is a single stratum
    JK.gb = D.gbase[D.strat, 1]
    // positions in the cumulative arrays (one leading zero per stratum)
    JK.cp1 = JK.gb + D.strat + D.rm1
    JK.cp0 = JK.gb + D.strat + D.rm0
    // records at risk at one or more event times, by their last event at risk
    act = (D.rm1 :> D.rm0)
    idx = selectindex(act)
    if (rows(idx) == 0) _error(2000, "no records at risk at the event times")
    _stpaf_jkpanel(idx, JK.gb, D.rm1, X, o, info, k, Xs)
    JK.o1 = o
    JK.in1 = info
    JK.k1 = k
    JK.X1 = Xs
    // ... and by the last event before they enter the risk set
    idx = selectindex(act :& (D.rm0 :> 0))
    _stpaf_jkpanel(idx, JK.gb, D.rm0, X, o, info, k, Xs)
    JK.o0 = o
    JK.in0 = info
    JK.k0 = k
    JK.X0 = Xs
    // failures (records with zero weight in the full sample never define event times)
    idx = selectindex((D.d :!= 0) :& (D.w :> 0))
    _stpaf_jkpanel(idx, JK.gb, D.rm1, X, o, info, k, Xs)
    JK.of = o
    JK.inf = info
    JK.kf = k
    JK.Xf = Xs
    return(JK)
}

// sums of the rows of V (sorted by event index) into the Ktot event times
real matrix _stpaf_bsum(real matrix V, real matrix info, real colvector k, real scalar K)
{
    real matrix out

    out = J(K, cols(V), 0)
    if (rows(info)) out[k, .] = panelsum(V, info)
    return(out)
}

// risk-set sums at every event time (sum over records with m0 < event <= m1) of the
// rows of V1 (records sorted as o1) and V0 (the same records sorted as o0)
real matrix _stpaf_jkrisk(struct stpaf_jk scalar JK, real matrix V1, real matrix V0)
{
    real matrix B, out
    real colvector v, r
    real scalar s, K, lo, hi, j

    B = _stpaf_bsum(V1, JK.in1, JK.k1, JK.Ktot)
    if (rows(JK.in0)) B = B - _stpaf_bsum(V0, JK.in0, JK.k0, JK.Ktot)
    // backward cumulation within strata
    out = B
    for (s = 1; s <= JK.S; s++) {
        K = JK.Kst[s]
        if (K == 0) continue
        lo = JK.gbase[s] + 1
        hi = JK.gbase[s] + K
        for (j = 1; j <= cols(B); j++) {
            v = B[|lo, j \ hi, j|]
            r = quadrunningsum(v[(K::1), 1])
            out[|lo, j \ hi, j|] = r[(K::1), 1]
        }
    }
    return(out)
}

// forward cumulation with a leading zero per stratum (layout of D.Lc)
real colvector _stpaf_jkcum0(struct stpaf_jk scalar JK, real colvector v)
{
    real colvector out
    real scalar s, K, pb

    out = J(JK.Ktot + JK.S, 1, 0)
    for (s = 1; s <= JK.S; s++) {
        K = JK.Kst[s]
        if (K == 0) continue
        pb = JK.gbase[s] + s
        out[|pb + 1 \ pb + K|] = quadrunningsum(v[|JK.gbase[s] + 1 \ JK.gbase[s] + K|])
    }
    return(out)
}

// weighted Breslow log partial likelihood (up to a constant), score, information, and
// the cumulative baseline hazard Lc, at b, with record weights w
void _stpaf_jkll(struct stpaf_jk scalar JK, real matrix X, real colvector w, real rowvector b,
    real scalar ll, real rowvector U, real matrix I, real colvector Lc)
{
    real colvector eta, a, a1, a0, wf, dk, S0, S0s, pos, c, ac, ii
    real matrix S1, Eb, V0

    eta = X*b' :- JK.cen
    a = w :* exp(eta)
    wf = w[JK.of, 1]
    dk = _stpaf_bsum(wf, JK.inf, JK.kf, JK.Ktot)
    a1 = a[JK.o1, 1]
    if (rows(JK.o0)) {
        a0 = a[JK.o0, 1]
        V0 = a0 :* JK.X0
    }
    else {
        a0 = J(0, 1, .)
        V0 = J(0, JK.p, .)
    }
    S0 = _stpaf_jkrisk(JK, a1, a0)
    S1 = _stpaf_jkrisk(JK, a1 :* JK.X1, V0)
    pos = (dk :> 0)
    // event times whose failures all belong to the left-out cluster: dk = 0, no term
    S0s = S0 + (1 :- pos)
    Lc = _stpaf_jkcum0(JK, dk :/ S0s)
    c = Lc[JK.cp1, 1] - Lc[JK.cp0, 1]
    ac = a :* c
    ii = selectindex(pos)
    ll = quadsum(wf :* eta[JK.of, 1]) - quadsum(dk[ii, 1] :* ln(S0[ii, 1]))
    U = quadcolsum(wf :* JK.Xf) - quadcross(ac, X)
    Eb = S1 :/ S0s
    I = quadcross(X, ac, X) - quadcross(Eb, dk, Eb)
}

// Newton-Raphson with step halving; converged when mreldif(b, previous b) <= 1e-12
// and U*I^-1*U' <= 1e-12 (tight and scale-free, so the solution does not depend on
// the starting values); it returns the number of iterations, Lc that of the solution,
// and rc = 0, or 3498 (singular information) or 430 (no convergence, likelihood not
// finite), in which case the caller refits stcox for this cluster
real rowvector _stpaf_jknewton(struct stpaf_jk scalar JK, real matrix X, real colvector w,
    real rowvector b0, real scalar it, real colvector Lc, real scalar rc)
{
    real rowvector b, bp, bn, U, Un, step
    real matrix I, In, Ii
    real colvector Lcn
    real scalar ll, lln, nr, h

    if (quadsum(w[JK.of, 1]) <= 0) _error(2000, "no failures left in a jackknife replicate")
    rc = 0
    b = b0
    bp = b0
    it = 0
    _stpaf_jkll(JK, X, w, b, ll, U, I, Lc)
    if (ll >= .) {
        rc = 430
        return(b)
    }
    while (1) {
        Ii = invsym(I)
        if (diag0cnt(Ii)) {
            rc = 3498
            return(b)
        }
        step = U*Ii
        nr = step*U'
        if (nr >= .) {
            rc = 430
            return(b)
        }
        if (it > 0) {
            if (mreldif(b, bp) <= 1e-12 & nr <= 1e-12) break
        }
        if (it >= 100) {
            rc = 430
            return(b)
        }
        bn = b + step
        _stpaf_jkll(JK, X, w, bn, lln, Un, In, Lcn)
        h = 0
        while ((lln >= . | lln < ll - 1e-12*(abs(ll) + 1)) & h < 50) {
            step = step/2
            bn = b + step
            _stpaf_jkll(JK, X, w, bn, lln, Un, In, Lcn)
            h++
        }
        if (lln >= .) {
            rc = 430
            return(b)
        }
        bp = b
        b = bn
        ll = lln
        U = Un
        I = In
        Lc = Lcn
        it++
    }
    return(b)
}

// replicate PAFs for every scenario and horizon (columns in the row order of TAB);
// Q1, Q0: positions in Lc of each record's interval truncated at each horizon
real rowvector _stpaf_jkpaf(struct stpaf_jk scalar JK, real matrix X, real colvector w,
    real colvector g, real rowvector b, real colvector Lc, real scalar nsc,
    real matrix Q1, real matrix Q0)
{
    pointer(real matrix) scalar pX
    real rowvector out
    real colvector a, ac, dL
    real scalar j, t, nt, E, E0

    nt = cols(Q1)
    out = J(1, nsc*nt, .)
    a = w :* g :* exp(X*b' :- JK.cen)
    for (j = 1; j <= nsc; j++) {
        pX = findexternal("__stpaf_Xc" + strofreal(j))
        if (pX == NULL) _error("stpaf internal error: counterfactual design matrix not found")
        ac = w :* g :* exp((*pX)*b' :- JK.cen)
        for (t = 1; t <= nt; t++) {
            dL = Lc[Q1[., t], 1] - Lc[Q0[., t], 1]
            E = quadsum(a :* dL)
            E0 = quadsum(ac :* dL)
            if (E <= 0 | E0 <= 0) {
                _error(2000, "no expected events in the selected subset/horizon of a jackknife replicate")
            }
            out[1, (j-1)*nt + t] = 1 - E0/E
        }
    }
    return(out)
}

// driver of the fast jackknife; writes the mean number of Newton iterations (over the
// clusters fitted by Newton) to the Stata scalar itname, or -1 (and leaves TAB unchanged)
// when the model has tied failure times and was not fitted with Breslow ties, so that
// Stata refits stcox instead. nfname: number of clusters whose Newton fit failed; their
// codes go to the Stata matrix flname and the replicate PAFs of the others to the external
// __stpaf_PF (TAB unchanged), for Stata to refit those clusters and call _stpaf_jackfill
void _stpaf_jackfast(real matrix X, real colvector t0, real colvector t1, real colvector d,
    real colvector g, real colvector w, real colvector clus, real colvector strat,
    real rowvector b, real scalar nsc, real rowvector tau, real scalar level,
    real scalar logci, string scalar tabname, real scalar warm, real scalar brestie,
    string scalar itname, string scalar nfname, string scalar flname)
{
    struct stpaf_dat scalar D
    struct stpaf_jk scalar JK
    pointer() scalar pPF
    real matrix PF, Q1, Q0, mt
    real colvector cnt, wj, Lc, fl
    real rowvector bj, b0
    real scalar G, nt, s, t, gg, it, itsum, nok, rcg

    D = _stpaf_prep(X, t0, t1, d, g, w, clus, strat, b)
    JK = _stpaf_jkset(X, D)
    if (!brestie) {
        // Efron/exact likelihoods differ from Breslow only at tied failure times
        cnt = _stpaf_bsum(J(rows(JK.of), 1, 1), JK.inf, JK.kf, JK.Ktot)
        if (max(cnt) > 1) {
            st_numscalar(itname, -1)
            st_numscalar(nfname, 0)
            return
        }
    }
    G = D.G
    nt = cols(tau)
    mt = J(JK.S, nt, 0)
    for (s = 1; s <= JK.S; s++) {
        if (JK.Kst[s] == 0) continue
        for (t = 1; t <= nt; t++) mt[s, t] = _stpaf_pos(D.tk, JK.gbase[s], JK.Kst[s], tau[t])
    }
    Q1 = J(JK.n, nt, .)
    Q0 = J(JK.n, nt, .)
    for (t = 1; t <= nt; t++) {
        Q1[., t] = JK.gb + JK.strat + rowmin((JK.m1, mt[JK.strat, t]))
        Q0[., t] = JK.gb + JK.strat + rowmin((JK.m0, mt[JK.strat, t]))
    }
    b0 = (warm ? b : J(1, cols(b), 0))
    PF = J(G, nsc*nt, .)
    itsum = 0
    nok = 0
    fl = J(0, 1, .)
    for (gg = 1; gg <= G; gg++) {
        wj = w :* (clus :!= gg)
        bj = _stpaf_jknewton(JK, X, wj, b0, it, Lc, rcg)
        if (rcg) {
            fl = fl \ gg
            continue
        }
        itsum = itsum + it
        nok++
        PF[gg, .] = _stpaf_jkpaf(JK, X, wj, g, bj, Lc, nsc, Q1, Q0)
    }
    st_numscalar(itname, (nok > 0 ? itsum/nok : 0))
    st_numscalar(nfname, rows(fl))
    if (rows(fl)) {
        st_matrix(flname, fl)
        pPF = findexternal("__stpaf_PF")
        if (pPF == NULL) pPF = crexternal("__stpaf_PF")
        *pPF = PF
        return
    }
    _stpaf_jkvar(PF, level, logci, tabname)
}

// largest deviation (after removing a constant) between Stata's xb and ours
real scalar _stpaf_chk(struct stpaf_dat scalar D, string scalar xbvar, string scalar sel)
{
    real colvector x, df
    x = st_data(., xbvar, sel)
    df = x - D.xb
    return(max(abs(df :- mean(df))))
}

void _stpaf_post(struct stpaf_dat scalar D, string scalar sc, string scalar st)
{
    st_matrix(sc, D.scorechk)
    st_matrix(st, (D.G, D.n, sum(D.dk), D.Ktot, D.S))
}

end
