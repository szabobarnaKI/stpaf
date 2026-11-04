*! stpaf 0.14.0  04nov2026
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
*! 0.6.0: Newton jackknife with the Efron likelihood for efron models with tied failure times
*!        (stcox refits only remain for exactm/exactp ties with tied failure times)
*! 0.7.0: jackknife confidence intervals use the t distribution with G-1 degrees of freedom;
*!        fix: jackknife refits after exactm/exactp used breslow ties (e(ties) is marginal/partial)
*! 0.8.0: multiple imputation: after mi estimate: stcox, stcox is refitted in each imputation,
*!        the PAF and its variance are estimated there and combined by Rubin's rules (on the
*!        scale of citrans()), with Barnard-Rubin degrees of freedom for the jackknife;
*!        r(scen#) returns the scenario labels
*! 0.9.0: influence function exact for efron ties with tied failure times: the coefficient
*!        term uses the Efron score residuals and Efron information (as stcox's esr and
*!        vce(cluster)); r(V_b) cluster-robust VCE of the coefficients, r(if_ties)
*! 0.10.0: after mi estimate, the per-imputation stcox fits saved by mi estimate, saving() (or
*!        the file of mi estimate using, or mifile()) are used instead of refitting stcox, when
*!        they match the stcox command, its estimation sample and the data (score at e(b));
*!        otherwise stcox is refitted (as before); options mirefit and mifile();
*!        r(mi_nsaved), r(mi_fits), r(mi_file)
*! 0.11.0: option varlevels(varname[=#]): level-specific PAFs of a factor variable (each level
*!        moved to the reference in turn); they add up to the joint PAF (all levels moved);
*!        r(levelvar), r(reflevel), r(levels), r(sumdiff); fix: design matrices and stcox refits
*!        of factor variables with 2+ non-base levels used a wrong base (stopped with 498);
*!        fix: counterfactual columns of continuous interactions (c.z#c.age) of a changed variable
*!        kept their observed values (fvrevar reused the observed temporaries; wrong PAFs)
*!        WARNING: results of versions before 0.11.0 are wrong when a variable changed by
*!        eliminate(), set() or cap() enters an interaction as a continuous variable
*! 0.12.0: option combined(factorlist): average sequential (Shapley; Eide and Gefeller) PAFs of
*!        K factors (2 <= K <= 8), linear combinations of the 2^K - 1 subset PAFs, which add up
*!        to the joint PAF; option subsets; r(K), r(factors), r(W_shap), r(subsets), r(sumdiff);
*!        r(V_paf): covariance of the PAFs of all rows of r(table) (every call, also after mi)
*! 0.13.0: memory and time: counterfactual design matrices are stored as column versions (only
*!        columns of terms that contain a changed variable, one version per distinct set of
*!        changes among the term's variables), built in as few passes over the data as possible
*!        (one for main effects); the influence-function covariance is accumulated over blocks
*!        of clusters (no clusters x rows matrix); vectorized influence function;
*!        r(xc_cols), r(xc_builds), r(if_chunks). Results equal 0.12.0 (points exactly,
*!        standard errors up to rounding)
*! 0.14.0: time: risk-set sums at the event times (observed and every counterfactual scenario)
*!        and the information and cluster scores computed with vectorized bucket sums instead of
*!        record loops; results equal 0.13.0 up to rounding

program define stpaf, rclass
    version 15.0
    // after mi estimate: stcox (also with mi estimate's post option, where e(cmd) is stcox):
    // per-imputation fits (saved or refitted) and Rubin's rules
    if ("`e(mi)'" == "mi" | "`e(cmd)'" == "mi estimate") & !strpos(`"`0'"', "misavedchk") {
        _stpaf_mi `0'
        return add
        exit
    }
    // misavedchk (internal, used by _stpaf_mi): stop with rc 7721 when e(b) does not maximize
    // the partial likelihood of the data in memory
    // ifchunk(#) (internal, for testing): clusters per block of the influence-function covariance
    syntax [if] [in] [, ELIMinate(varlist numeric) SET(string) CAP(string) ///
        HORizon(numlist >0 sort) Level(cilevel) CITrans(string) EACH VCE(string) NODOTS ///
        CLuster(varname) REFIT noWARM VARLEVels(string) COMBined(string) SUBsets MISAVEDCHK ///
        IFCHunk(integer 0) ]

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
    // ties of the model: 0 breslow, 1 efron, 2 exactm/exactp (stored by stcox as e(ties)
    // marginal/partial)
    local tiemode = cond(inlist("`e(ties)'", "breslow", ""), 0, cond("`e(ties)'" == "efron", 1, 2))

    // ------------------------------------------------- counterfactual scenario
    local nspec 0
    local lvar ""
    local K 0
    if `"`combined'"' != "" {
        // average sequential (Shapley) PAFs (0.12.0). Factors: each variable of an element
        // eliminate(), set() or cap() is one factor; the elements of a group in parentheses form
        // one factor. Internal scenarios: the 2^K - 1 non-empty subsets of the factors, the
        // joint scenario (all factors) first, then by size and lexicographic order. The Shapley
        // PAF of factor k is sum_S w_kS PAF(S), the average over the K! orders of changing the
        // factors of the extra PAF from also changing factor k; they add up to the joint PAF.
        if "`eliminate'`set'`cap'`each'" != "" | `"`varlevels'"' != "" {
            di as err "option {bf:combined()} may not be combined with {bf:eliminate()}, {bf:set()}, {bf:cap()}," ///
                " {bf:each} or {bf:varlevels()}"
            exit 198
        }
        // top-level elements (separated by blanks outside parentheses)
        mata: _stpaf_toplocals(st_local("combined"), "ctk")
        if "`ctkbad'" == "1" {
            di as err `"combined(): unbalanced parentheses in {bf:`combined'}"'
            exit 198
        }
        local ngrp 0
        local nspc 0
        forvalues i = 1/`ctkn' {
            local tk `"`ctk`i''"'
            if substr(`"`tk'"', 1, 1) == "(" {
                mata: st_local("isg", strofreal(_stpaf_isgroup(st_local("tk"))))
                if !`isg' {
                    di as err `"combined(): invalid element {bf:`tk'}"'
                    exit 198
                }
                local inner = substr(`"`tk'"', 2, strlen(`"`tk'"') - 2)
                local ++ngrp
                mata: _stpaf_toplocals(st_local("inner"), "gtk")
                if "`gtkbad'" == "1" | `gtkn' == 0 {
                    di as err `"combined(): invalid group {bf:`tk'}"'
                    exit 198
                }
                forvalues jj = 1/`gtkn' {
                    local ++nspc
                    local cspec`nspc' `"`gtk`jj''"'
                    local cgrp`nspc' `ngrp'
                }
            }
            else {
                local ++nspc
                local cspec`nspc' `"`tk'"'
                local cgrp`nspc' 0
            }
        }
        // elements -> changes (sv, sc, sl) and factors (fit = changes, flab = label)
        local allv ""
        forvalues j = 1/`nspc' {
            local sp `"`cspec`j''"'
            if !regexm(`"`sp'"', "^([A-Za-z]+)[(](.*)[)]$") {
                di as err `"combined(): invalid element {bf:`sp'}; use {bf:eliminate()}, {bf:set()} or {bf:cap()},"' ///
                    " or several of them in parentheses (one level) for a factor"
                exit 198
            }
            local onm = lower(regexs(1))
            local oct = strtrim(regexs(2))
            local kind ""
            if "`onm'" == "set" local kind "set"
            else if "`onm'" == "cap" local kind "cap"
            else if strlen("`onm'") >= 4 & "`onm'" == substr("eliminate", 1, strlen("`onm'")) local kind "elim"
            if "`kind'" == "" | `"`oct'"' == "" {
                di as err `"combined(): invalid element {bf:`sp'}; use {bf:eliminate()}, {bf:set()} or {bf:cap()}"'
                exit 198
            }
            if "`kind'" == "elim" {
                unab vl : `oct', name(combined())
                local items ""
                foreach v of local vl {
                    confirm numeric variable `v'
                    local items `items' `v'=0
                }
            }
            else {
                local oct = stritrim(`"`oct'"')
                local oct : subinstr local oct " =" "=", all
                local oct : subinstr local oct "= " "=", all
                local items ""
                foreach tok of local oct {
                    if !regexm("`tok'", "^([A-Za-z_][A-Za-z0-9_]*)=(.+)$") {
                        di as err "combined(): invalid specification {bf:`tok'} in {bf:`sp'}; use varname=#"
                        exit 198
                    }
                    local v = regexs(1)
                    local val = regexs(2)
                    unab v : `v', max(1) name(combined())
                    confirm numeric variable `v'
                    confirm number `val'
                    local items `items' `v'=`val'
                }
            }
            foreach it of local items {
                local eqp = strpos("`it'", "=")
                local v = substr("`it'", 1, `eqp' - 1)
                local val = substr("`it'", `eqp' + 1, .)
                local isdup : list v in allv
                if `isdup' {
                    di as err "combined(): variable {bf:`v'} appears more than once"
                    exit 198
                }
                local allv `allv' `v'
                local ++nspec
                local sv`nspec' "`v'"
                if "`kind'" == "cap" {
                    local sc`nspec' "qui replace `v' = min(`v', `val')"
                    local sl`nspec' "`v'<=`val'"
                }
                else {
                    local sc`nspec' "qui replace `v' = `val'"
                    local sl`nspec' "`v'=`val'"
                }
                // factor: the group's (numbered at its first change) or a new one per variable
                if `cgrp`j'' > 0 {
                    if "`gf`cgrp`j'''" == "" {
                        local ++K
                        local gf`cgrp`j'' `K'
                    }
                    local f `gf`cgrp`j'''
                }
                else {
                    local ++K
                    local f `K'
                }
                local fit`f' `fit`f'' `nspec'
                if "`flab`f''" == "" local flab`f' "`v'"
                else local flab`f' "`flab`f''+`v'"
            }
        }
        if `K' < 2 {
            di as err "combined() requires at least two factors; for one factor use {bf:eliminate()}, {bf:set()} or {bf:cap()}"
            exit 198
        }
        if `K' > 8 {
            di as err "combined() allows at most 8 factors (`K' given; 2^K - 1 subset scenarios are computed)"
            exit 198
        }
        tempname SUBM WSH
        mata: _stpaf_shapmat(`K', "`SUBM'", "`WSH'")
        local nsc = 2^`K' - 1
        forvalues s = 1/`nsc' {
            local scl`s' ""
            local lab`s' ""
            forvalues f = 1/`K' {
                if `SUBM'[`s', `f'] == 1 {
                    local scl`s' `scl`s'' `fit`f''
                    if "`lab`s''" == "" local lab`s' "`flab`f''"
                    else local lab`s' "`lab`s'',`flab`f''"
                }
            }
        }
        local lab1 "joint"
        local jl ""
        forvalues i = 1/`nspec' {
            local jl `jl' `sl`i''
        }
        local flist ""
        forvalues f = 1/`K' {
            local flist `flist' `flab`f''
        }
    }
    else if `"`varlevels'"' != "" {
        // level-specific PAFs of a factor variable (0.11.0): scenario 1 moves every record to
        // the reference level (joint), scenario 1+k only the records at the k-th other level.
        // Each record has one level, so E - E_cf of the joint scenario is the sum over the
        // levels: the level PAFs add up to the joint PAF.
        if "`eliminate'`set'`cap'`each'" != "" {
            di as err "option {bf:varlevels()} may not be combined with {bf:eliminate()}, {bf:set()}, {bf:cap()} or {bf:each}"
            exit 198
        }
        local vlspec = subinstr(`"`varlevels'"', " ", "", .)
        if !regexm(`"`vlspec'"', "^([A-Za-z_][A-Za-z0-9_]*)(=([0-9]+))?$") {
            di as err `"varlevels(): invalid specification {bf:`varlevels'};"' ///
                " use {it:varname} or {it:varname}={it:#} with # a level of the factor variable"
            exit 198
        }
        local lvar = regexs(1)
        // (regexs() of an optional group that did not match prints an error message)
        local lref ""
        local eqp = strpos(`"`vlspec'"', "=")
        if `eqp' local lref = substr(`"`vlspec'"', `eqp' + 1, .)
        unab lvar : `lvar', max(1) name(varlevels())
        confirm numeric variable `lvar'
        // levels of the variable in the factor-variable terms of e(b) (also omitted ones: their
        // records exist and the interactions change), and the base level
        local cnall : colnames e(b)
        local lvls ""
        local lbase ""
        foreach t of local cnall {
            local parts = subinstr("`t'", "#", " ", .)
            foreach pt of local parts {
                if regexm("`pt'", "^([0-9]+)([a-z]*)\.(.+)$") {
                    local plv = regexs(1)
                    local pfl = regexs(2)
                    local pvn = regexs(3)
                    if "`pvn'" == "`lvar'" {
                        local plv = real("`plv'")
                        local lvls `lvls' `plv'
                        if strpos("`pfl'", "b") & !strpos("`pfl'", "bn") local lbase `lbase' `plv'
                    }
                }
            }
        }
        if "`lvls'" == "" {
            di as err "varlevels(): {bf:`lvar'} does not enter the model as a factor variable ({bf:i.`lvar'})"
            exit 198
        }
        local lvls : list uniq lvls
        numlist "`lvls'", sort
        local lvls "`r(numlist)'"
        local lbase : list uniq lbase
        if "`lref'" == "" {
            if `: word count `lbase'' != 1 {
                di as err "varlevels(): the base level of {bf:`lvar'} cannot be determined from e(b);" ///
                    " specify the reference level as {bf:varlevels(`lvar'=}{it:#}{bf:)}"
                exit 198
            }
            local lref `lbase'
        }
        else {
            local lref = real("`lref'")
            local inl : list lref in lvls
            if !`inl' {
                di as err "varlevels(): `lref' is not a level of {bf:`lvar'} in the model (levels: `lvls')"
                exit 198
            }
        }
        local lvnr : list lvls - lref
        if "`lvnr'" == "" {
            di as err "varlevels(): {bf:`lvar'} has no level other than the reference in the model"
            exit 198
        }
        local nspec 1
        local sv1 "`lvar'"
        local sc1 "qui replace `lvar' = `lref'"
        local sl1 "`lvar'=`lref'"
        local nsc 1
        local scl1 "1"
        local lab1 "`lvar'=`lref'"
        foreach k of local lvnr {
            local ++nspec
            local sv`nspec' "`lvar'"
            local sc`nspec' "qui replace `lvar' = `lref' if `lvar' == `k'"
            local sl`nspec' "`k'.`lvar'"
            local ++nsc
            local scl`nsc' "`nspec'"
            local lab`nsc' "`k'.`lvar'"
        }
        local jl "`lvar'=`lref', each level"
    }
    else {
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
        di as err "specify at least one of eliminate(), set(), cap(), varlevels() or combined()"
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
    }

    // ------------------------------------------------------- samples, ids, strata
    tempvar esamp clus strat xbchk wvar
    qui gen byte `esamp' = e(sample)
    marksample touse, novarlist
    qui replace `touse' = 0 if `esamp' != 1
    qui count if `touse'
    if r(N) == 0 error 2000

    // clusters: option cluster() > e(clustvar) > stset id() > single records.
    // The variance is built here from the partial-likelihood information and the cluster
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
    // design matrix from ALL column names of e(b) (with base and omitted terms such as 0b.x),
    // then the non-omitted columns: given only 1.x 2.x 3.x, fvrevar takes the lowest level
    // listed (1) as the base and returns a column of zeros for 1.x (fixed in 0.11.0)
    fvrevar `cn'
    local xall "`r(varlist)'"
    local xvars ""
    foreach j of local keepj {
        local xvars `xvars' `: word `j' of `xall''
    }
    qui predict double `xbchk' if `esamp', xb

    // ------------------------- counterfactual design matrices as column versions (0.13.0)
    // A column of the design matrix (a term of e(b)) depends only on the changes of the
    // variables in the term (read from its name), so in scenario s it equals the version built
    // with the changes of s that act on these variables. Only such versions are stored, in the
    // Mata matrix __stpaf_XV (records x versions); __stpaf_XM (scenarios x columns) gives the
    // version of each column in each scenario (0 = observed column) and __stpaf_XC the column
    // of each version. A pass over the data with the changes of scenario s applied (preserve,
    // changes, fvrevar, restore) gives the versions of all columns for s at once; the passes
    // are chosen by _stpaf_xcplan (all versions with as few passes as possible, the joint
    // scenario first). In each pass the columns that should not change are compared with the
    // observed ones (internal check). The observed design matrix is kept in __stpaf_X0.
    // With efron ties and tied failure times (unit weights) the coefficient part of the
    // influence function uses the Efron score residuals and information (IFT = 1).
    // VP: covariance of the PAFs of all rows (influence function; replaced by the jackknife)
    tempname TAB SC ST CHK VB IFT VP XCK NCH
    local svall ""
    forvalues i = 1/`nspec' {
        unab vu : `sv`i'', max(1)
        local svall `svall' `vu'
    }
    qui count if `esamp'
    local nes = r(N)
    foreach x in XV XM XC X0 {
        capture mata: mata drop __stpaf_`x'
    }
    local xcnb 0
    local xcnu 0
    scalar `NCH' = .
    capture noisily {
        mata: _stpaf_xcplan("`fvl'", "`svall'", `nsc', `nes', "xc")
        mata: __stpaf_X0 = st_data(., tokens("`xvars'"), "`esamp'")
        forvalues b = 1/`xcnb' {
            local sb `xcs`b''
            preserve
            // drop the temporaries that fvrevar made for the observed design matrix: fvrevar
            // reuses an existing temporary for a term without noticing that the data changed
            // (c.z#c.age kept its observed values after z was changed; fixed in 0.11.0)
            foreach v of local xall {
                if substr("`v'", 1, 2) == "__" capture drop `v'
            }
            foreach i of local scl`sb' {
                capture recast double `sv`i''
                `sc`i''
            }
            fvrevar `cn'
            local xca "`r(varlist)'"
            local vn ""
            foreach jp of local xcc`b' {
                local jj : word `jp' of `keepj'
                local vn `vn' `: word `jj' of `xca''
            }
            local kn ""
            foreach jp of local xck`b' {
                local jj : word `jp' of `keepj'
                local kn `kn' `: word `jj' of `xca''
            }
            mata: _stpaf_xcput("`vn'", "`xcv`b''", "`kn'", "`xck`b''", "`esamp'", "`XCK'")
            restore
            if scalar(`XCK') > 0 {
                local bad : word `=scalar(`XCK')' of `fvl'
                di as err "stpaf internal error: term {bf:`bad'} changes in a counterfactual scenario" ///
                    " although it contains none of the changed variables"
                exit 498
            }
        }
        mata: _stpaf_run(__stpaf_X0, ///
            st_data(., "_t0", "`esamp'"), st_data(., "_t", "`esamp'"), ///
            st_data(., "_d", "`esamp'"), st_data(., "`touse'", "`esamp'"), ///
            st_data(., "`wvar'", "`esamp'"), ///
            st_data(., "`clus'", "`esamp'"), st_data(., "`strat'", "`esamp'"), ///
            st_matrix("`bk'"), `nsc', (`tauvec'), `level', `logci', `fwrep', ///
            "`TAB'", "`SC'", "`ST'", "`CHK'", "`xbchk'", "`esamp'", ///
            `tiemode', "`VB'", "`IFT'", "`VP'", `ifchunk', "`NCH'")
    }
    local rc = _rc
    if !`rc' & scalar(`CHK') > 1e-6 {
        di as err "stpaf could not reproduce the linear predictor of stcox (max abs diff = " ///
            %9.3g scalar(`CHK') ")"
        di as err "the model probably contains terms that cannot be rebuilt from e(b)"
        local rc 498
    }
    local ifties "breslow"
    // largest standardized partial-likelihood score at e(b) (Efron when IFT = 1, else Breslow)
    local maxsc = 0
    if !`rc' {
        if scalar(`IFT') == 1 local ifties "efron"
        forvalues j = 1/`p' {
            if `SC'[1, `j'] > `maxsc' local maxsc = `SC'[1, `j']
        }
        // a fit saved by mi estimate is used only if e(b) maximizes the partial likelihood
        // of the data in memory (otherwise _stpaf_mi refits stcox in this imputation)
        if "`misavedchk'" != "" & !(`maxsc' <= 1e-4) local rc 7721
    }

    // ----------------------------------------------------- exact cluster jackknife
    // Each cluster is left out in turn and the model is re-estimated on the reduced
    // sample (same terms, strata and stset weights); the PAFs of the reduced samples
    // give the jackknife variance of log(1-PAF) and of PAF, which replaces the
    // influence-function one.
    // Default (breslow or efron ties, or exactm/exactp without tied failure times):
    // Newton-Raphson for the Breslow or (efron with tied failure times) Efron partial
    // likelihood in Mata, the left-out cluster given zero weight, started from e(b) (nowarm:
    // from 0) and iterated until mreldif(b) <= 1e-12 and U'I^-1U <= 1e-12, i.e. to convergence
    // independent of the starting values.
    // Clusters whose Newton fit fails (singular information, no convergence) are refitted
    // by stcox. Option refit, exactm/exactp ties with tied failure times, and weighted
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
            local dowarm = ("`warm'" == "")
            if `tiemode' & "`wtype'" != "" {
                // defensive: stcox 19 allows no weights with efron (cert 12d); if a weighted
                // efron/exact fit occurs, fweights > 1 imply tied failures that a count of
                // failure records cannot see, so refit with the model's own ties
                di as txt "note: weighted model with {bf:`e(ties)'} ties; the jackknife refits {bf:stcox}"
                if "`warm'" != "" di as txt "      (option nowarm has no effect)"
            }
            else {
                `capn' mata: _stpaf_jackfast(__stpaf_X0, ///
                    st_data(., "_t0", "`esamp'"), st_data(., "_t", "`esamp'"), ///
                    st_data(., "_d", "`esamp'"), st_data(., "`touse'", "`esamp'"), ///
                    st_data(., "`wvar'", "`esamp'"), ///
                    st_data(., "`clus'", "`esamp'"), st_data(., "`strat'", "`esamp'"), ///
                    st_matrix("`bk'"), `nsc', (`tauvec'), `level', `logci', "`TAB'", ///
                    `dowarm', `tiemode', "`JKI'", "`JKNF'", "`JKF'", "`VP'")
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
                        di as txt "note: tied failure times with {bf:`e(ties)'} ties (no Newton jackknife for" ///
                            " the exact likelihoods); the jackknife refits {bf:stcox}"
                        if "`warm'" != "" di as txt "      (option nowarm has no effect)"
                    }
                }
            }
        }
        // stcox refits: all clusters (full) or those whose Newton fit failed (JKF)
        local nref = cond(`full', `G0', `jknref')
        if !`rc' & `nref' > 0 {
            // stcox stores e(ties) = "marginal" after exactm and "partial" after exactp
            local ties "`e(ties)'"
            if "`ties'" == "marginal" local ties "exactm"
            else if "`ties'" == "partial" local ties "exactp"
            else if !inlist("`ties'", "breslow", "efron", "exactm", "exactp") local ties "breslow"
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
                capture quietly stcox `cn' if `esamp' & `clus' != `g', `stopt' `ties'
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
                    `capn' mata: _stpaf_jack(__stpaf_X0, ///
                        st_data(., "_t0", "`esamp'"), st_data(., "_t", "`esamp'"), ///
                        st_data(., "_d", "`esamp'"), st_data(., "`touse'", "`esamp'"), ///
                        st_data(., "`wvar'", "`esamp'"), ///
                        st_data(., "`clus'", "`esamp'"), st_data(., "`strat'", "`esamp'"), ///
                        st_matrix("`BJ'"), `nsc', (`tauvec'), `level', `logci', "`TAB'", "`VP'")
                }
                else {
                    `capn' mata: _stpaf_jackfill(__stpaf_X0, ///
                        st_data(., "_t0", "`esamp'"), st_data(., "_t", "`esamp'"), ///
                        st_data(., "_d", "`esamp'"), st_data(., "`touse'", "`esamp'"), ///
                        st_data(., "`wvar'", "`esamp'"), ///
                        st_data(., "`clus'", "`esamp'"), st_data(., "`strat'", "`esamp'"), ///
                        st_matrix("`BJ'"), st_matrix("`JKF'"), `nsc', (`tauvec'), `level', `logci', ///
                        "`TAB'", "`VP'")
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
    foreach x in XV XM XC X0 {
        capture mata: mata drop __stpaf_`x'
    }
    if `rc' exit `rc'

    // ------------------------------------------- average sequential (Shapley) PAFs
    // fixed linear combinations of the subset PAFs (point estimates, influence functions and
    // jackknife replicates alike): Shapley rows inserted after the joint scenario; their
    // standard errors from the covariance VP, CIs on the PAF scale
    local nscin = `nsc'
    if `K' > 0 {
        local dfq = cond("`vce'" == "jackknife", `ST'[1, 1] - 1, .)
        mata: _stpaf_shapfin("`TAB'", "`VP'", "`WSH'", `nt', `level', `dfq')
        forvalues s = 1/`nscin' {
            local ilab`s' "`lab`s''"
        }
        local lab1 "joint"
        forvalues k = 1/`K' {
            local lab`=`k'+1' "Shapley: `flab`k''"
        }
        forvalues s = 2/`nscin' {
            local lab`=`K'+`s'' "`ilab`s''"
        }
        local nsc = `nscin' + `K'
    }

    // ---------------------------------------------------------------- display
    local G    = `ST'[1, 1]
    local nev  = `ST'[1, 3]
    local nK   = `ST'[1, 4]
    local nS   = `ST'[1, 5]

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
            local jkdesc = "`=`G'-`jknref'' Newton fits" + cond("`e(ties)'" == "efron", " (Efron)", "") + ", " + ///
                cond("`warm'" == "", "warm start", "cold start") + ///
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
        di as txt "Std. err. from the influence function" ///
            cond("`ifties'" == "efron", " (Efron score residuals)", "") ", clustered on " ///
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
        if `K' > 0 {
            // subset rows only with option subsets; group headings
            if `sj' > `K' + 1 & "`subsets'" == "" continue
            if mod(`i' - 1, `nt') == 0 {
                if `sj' == 2 di as txt "{ralign 12:Shapley} {c |}"
                else if `sj' == `K' + 2 di as txt "{ralign 12:Subsets} {c |}"
            }
        }
        local dl : subinstr local lab`sj' "Shapley: " ""
        di as txt "{ralign 12:" abbrev("`dl'", 12) "} {c |}" as res _col(17) %7.4g `h' ///
            _col(26) %8.4f `TAB'[`i', 3] _col(36) %8.4f `TAB'[`i', 4] ///
            _col(48) "[" %7.4f `TAB'[`i', 5] "," %7.4f `TAB'[`i', 6] "]" ///
            _col(70) %8.1f `TAB'[`i', 7] _col(80) %8.1f `TAB'[`i', 8]
    }
    di as txt "{hline 13}{c BT}{hline 70}"
    // level PAFs and Shapley PAFs: largest |sum of the rows - joint PAF| over the horizons
    // (exact up to rounding)
    tempname SDIF SSUM
    scalar `SDIF' = .
    local s1 0
    if "`lvar'" != "" {
        local s1 2
        local s2 `nsc'
    }
    if `K' > 0 {
        local s1 2
        local s2 = `K' + 1
    }
    if `s1' {
        scalar `SDIF' = 0
        forvalues h = 1/`nt' {
            scalar `SSUM' = 0
            forvalues j = `s1'/`s2' {
                scalar `SSUM' = `SSUM' + `TAB'[(`j'-1)*`nt' + `h', 3]
            }
            scalar `SDIF' = max(`SDIF', abs(`SSUM' - `TAB'[`h', 3]))
        }
    }
    if "`lvar'" != "" {
        di as txt "Level PAFs ({it:k}.`lvar'): records at level {it:k} moved to `lvar' = `lref', other records"
        di as txt "unchanged; they add up to the joint PAF (`lvar'=`lref')."
    }
    if `K' > 0 {
        di as txt "Shapley PAFs: the extra PAF from also changing the factor, averaged over all `K'! orders"
        di as txt "of changing the factors (Eide and Gefeller 1995); they add up to the joint PAF. Their"
        di as txt "CIs are on the PAF scale" _c
        if "`vce'" == "jackknife" di as txt ", with t(" `G'-1 ") quantiles." _c
        else di as txt "." _c
        if "`subsets'" == "" {
            di as txt " The " 2^`K'-2 " other subset scenarios are in r(table)"
            di as txt "(option {bf:subsets} displays them)."
        }
        else di
    }
    local ciwho = cond(`K' > 0, "Joint and subset rows: CI", "CI")
    if "`citrans'" == "log" {
        di as txt "`ciwho' is computed on the log(1-PAF) scale and back-transformed" _c
        if "`vce'" == "jackknife" di as txt ", with t(" `G'-1 ") quantiles."
        else di as txt "."
    }
    else if "`vce'" == "jackknife" {
        di as txt "`ciwho' uses t(" `G'-1 ") quantiles."
    }
    if "`vce'" == "if" & `G' < 100 & !`jkfail' & !`fwrep' {
        di as txt "Note: with fewer than 100 clusters the influence-function CI tends to be too narrow;"
        di as txt "      {bf:vce(jackknife)}, the default in this case, is recommended."
    }
    // left: exactm/exactp ties with tied failure times (or a weighted efron fit)
    if `maxsc' > 1e-4 & "`vce'" == "if" {
        di as txt "Note: standardized Breslow partial-likelihood score at e(b) = " %9.2g `maxsc' ///
            "; the model was probably fitted with exactm or exactp ties."
        di as txt "      Breslow-type expressions are used for the coefficient part of the variance;" ///
            " {bf:vce(jackknife)} or refitting with efron or breslow ties is exact."
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
    // storage diagnostics: counterfactual column versions stored, passes over the data,
    // blocks of clusters of the influence-function covariance
    return scalar xc_cols   = `xcnu'
    return scalar xc_builds = `xcnb'
    return scalar if_chunks = scalar(`NCH')
    // covariance of the PAFs of the rows of r(table) (row names s#h#: scenario, horizon)
    local vpn ""
    forvalues i = 1/`nr' {
        local vpn `vpn' s`=`TAB'[`i', 1]'h`=mod(`i' - 1, `nt') + 1'
    }
    matrix rownames `VP' = `vpn'
    matrix colnames `VP' = `vpn'
    return matrix V_paf = `VP'
    matrix colnames `TAB' = scenario horizon paf se lb ub E E_cf se_logR events
    return matrix table = `TAB'
    // cluster-robust VCE of the (non-omitted) coefficients used by the influence function
    capture confirm matrix `VB'
    if !_rc {
        matrix colnames `VB' = `fvl'
        matrix rownames `VB' = `fvl'
        return matrix V_b = `VB'
    }
    return local if_ties "`ifties'"
    return local cfactual "`jl'"
    return local citrans "`citrans'"
    return local vce "`vce'"
    return local clustvar "`clv'"
    forvalues j = 1/`nsc' {
        return local scen`j' "`lab`j''"
    }
    if "`lvar'" != "" {
        return local levelvar "`lvar'"
        return scalar reflevel = `lref'
        return local levels "`lvnr'"
        return scalar sumdiff = `SDIF'
    }
    if `K' > 0 {
        // columns of W_shap and rows of subsets: the subset scenarios in r(table) order
        // (scenario 1 = all factors, then scenarios K+2, ...); factors f1..fK in r(factors) order
        local wn "scen1"
        forvalues s = 2/`nscin' {
            local wn `wn' scen`=`K'+`s''
        }
        local fn ""
        forvalues f = 1/`K' {
            local fn `fn' f`f'
        }
        matrix rownames `WSH' = `fn'
        matrix colnames `WSH' = `wn'
        matrix rownames `SUBM' = `wn'
        matrix colnames `SUBM' = `fn'
        return matrix W_shap = `WSH'
        return matrix subsets = `SUBM'
        return scalar K = `K'
        return local factors "`flist'"
        return scalar sumdiff = `SDIF'
    }
    if "`vce'" == "jackknife" {
        return local jk_method "`jkmeth'"
        return scalar df = `G' - 1
        if "`jkmeth'" != "refit" {
            return scalar jk_iter = `jkit'
            return scalar jk_nrefit = `jknref'
        }
    }
    return local cmd "stpaf"
end

// =============================================================================
// Multiple imputation (0.8.0; saved fits 0.10.0)
// After mi estimate: stcox, stpaf is run in each imputation used by mi estimate (e(m_est_mi))
// with the same options, and the results are combined by Rubin's rules on the CI scale:
// log(1-PAF) (citrans(log), the default) or PAF (citrans(none)); Shapley rows of combined()
// always on the PAF scale (0.12.0). Degrees of freedom:
// large-sample (Rubin 1987) with the influence function, Barnard-Rubin (1999) with
// complete-data df G-1 with the jackknife (as mi estimate does for commands that store e(df_r)).
// The completed-data stcox fits (0.10.0): from the file saved by mi estimate, saving() (or the
// file of mi estimate using, or mifile()), when a saved fit is of the same stcox command
// (e(cmdline)), its estimation sample (rebuilt from the command) has e(N) observations, and its
// e(b) maximizes the partial likelihood of the imputed data in memory (largest standardized
// score <= 1e-4, checked by stpaf, misavedchk); otherwise stcox is refitted in that imputation
// (the completed-data command is taken from e(cmdline_mi)). The mean of the coefficients used
// must reproduce e(b_mi) (data unchanged since mi estimate).
// r(V_paf) (0.12.0): mean within-imputation covariance + (1 + 1/M) between-imputation
// covariance of the PAFs, on the PAF scale.
// =============================================================================
program define _stpaf_mi, rclass
    version 15.0
    syntax [if] [in] [, Level(cilevel) CITrans(string) NODOTS MIREFIT MIFILE(string) SUBsets *]
    if "`citrans'" == "" local citrans "log"
    if !inlist("`citrans'", "log", "none") {
        di as err "citrans() must be {bf:log} or {bf:none}"
        exit 198
    }
    local logci = ("`citrans'" == "log")
    if "`e(cmd_mi)'" != "stcox" {
        di as err "stpaf after {bf:mi estimate} requires {bf:mi estimate: stcox}"
        exit 301
    }
    // completed-data command: the text after the first colon outside parentheses and quotes
    local micl `"`e(cmdline_mi)'"'
    mata: st_local("stcmd", _stpaf_micmd(st_local("micl")))
    // file of the per-imputation fits: mifile(), else saving() or using of mi estimate
    local mif `"`mifile'"'
    local mifsrc "mifile()"
    if `"`mif'"' == "" {
        mata: st_local("mif", _stpaf_mifile(st_local("micl")))
        local mifsrc "mi estimate"
    }
    if `"`mif'"' != "" mata: st_local("mif", _stpaf_sterfn(st_local("mif")))
    local havef 0
    local fmiss 0
    if `"`mif'"' != "" {
        capture confirm file `"`mif'"'
        if !_rc local havef 1
        else if `"`mifile'"' != "" {
            di as err `"file {bf:`mif'} not found"'
            exit 601
        }
        else local fmiss 1
    }
    local mlist "`e(m_est_mi)'"
    if "`mlist'" == "" local mlist "`e(m_mi)'"
    local M : word count `mlist'
    if `M' < 2 {
        di as err "at least 2 imputations are required"
        exit 2001
    }
    local esv = (e(esampvary_mi) == 1)
    local mfirst : word 1 of `mlist'
    tempname BMI BSUM ALL TABm H OUT PF SE VPS VP WSH SUBM
    matrix `BMI' = e(b_mi)
    quietly mi query
    if "`r(style)'" == "" {
        di as err "data are not {bf:mi set}"
        exit 119
    }
    // the default jackknife (fewer than 100 clusters) may fall back to the influence function
    local vceopt = (strpos(`"`options'"', "vce(") > 0)

    _estimates hold `H', restore
    // mi estimate using: the stcox command is taken from the saved fits
    if `"`stcmd'"' == "" & `havef' {
        capture quietly estimates use `"`mif'"', number(1)
        if !_rc & "`e(cmd2)'" == "stcox" local stcmd `"`e(cmdline)'"'
        if `"`stcmd'"' != "" mata: st_local("stcmd", _stpaf_cmdnorm(st_local("stcmd")))
    }
    gettoken w1 : stcmd
    if "`w1'" != "stcox" {
        di as err "the {bf:stcox} command could not be found in the {bf:mi estimate} command line:"
        di as err `"  `micl'"'
        if `fmiss' di as err `"(file {bf:`mif'} of the saved fits not found; see option {bf:mifile()})"'
        exit 198
    }
    // saved fits: number k of the file belongs to the k-th imputation of e(m_est_mi) when the
    // file holds exactly M fits, otherwise to imputation k
    local usef = (`havef' & "`mirefit'" == "")
    local posj 0
    if `usef' {
        // (the file of mi estimate, saving() holds the M fits followed by one record without
        // e(b), seen in Stata 19: a fit is a loaded result with e(b) and e(cmd2) stcox)
        capture quietly estimates use `"`mif'"', number(`M')
        local ok1 = !_rc
        if `ok1' {
            capture confirm matrix e(b)
            local ok1 = (!_rc & "`e(cmd2)'" == "stcox")
        }
        capture quietly estimates use `"`mif'"', number(`=`M'+1')
        local ok2 = !_rc
        if `ok2' {
            capture confirm matrix e(b)
            local ok2 = (!_rc & "`e(cmd2)'" == "stcox")
        }
        local posj = (`ok1' & !`ok2')
        mata: st_local("stnorm", _stpaf_cmdnorm(st_local("stcmd")))
    }
    tempvar es ww
    local nsaved 0
    local nld 0
    local ncm 0
    local nsm 0
    local nbs 0
    local ner 0
    local cmdx ""

    preserve
    if "`nodots'" == "" _dots 0, title(Imputations) reps(`M')
    local j 0
    local nfail 0
    local nref 0
    foreach m of local mlist {
        local ++j
        restore, preserve
        quietly mi extract `m', clear
        local saved 0
        if `usef' {
          // fit number in the file: the k-th of e(m_est_mi) when the file holds exactly M fits,
          // else number m; the other numbering is tried when the first fails a check (wrong
          // data -> the score check rejects it)
          local cands = cond(`posj', "`j' `m'", "`m' `j'")
          local cands : list uniq cands
          local why ""
          foreach pos of local cands {
            if `saved' continue, break
            capture quietly estimates use `"`mif'"', number(`pos')
            local rcl = _rc
            if !`rcl' {
                capture confirm matrix e(b)
                local rcl = _rc
            }
            if `rcl' local why "ld"
            else {
                mata: st_local("ecl", _stpaf_cmdnorm(st_global("e(cmdline)")))
                if "`e(cmd2)'" != "stcox" | `"`ecl'"' != `"`stnorm'"' {
                    local why "cm"
                    if `"`cmdx'"' == "" local cmdx `"`e(cmdline)'"'
                }
                else {
                    // estimation sample of the saved fit, rebuilt from the stcox command
                    capture drop `es'
                    quietly gen byte `es' = 0
                    capture _stpaf_misamp `es' `stcmd'
                    local rc = _rc
                    if !`rc' {
                        quietly count if `es'
                        local nes = r(N)
                        if "`e(wtype)'" == "fweight" {
                            local fwexp `"`e(wexp)'"'
                            if `"`fwexp'"' == "" {
                                local fwv : char _dta[st_wv]
                                local fwexp "= `fwv'"
                            }
                            capture drop `ww'
                            capture quietly gen double `ww' `fwexp' if `es'
                            if _rc local nes = .
                            else {
                                quietly summarize `ww' if `es', meanonly
                                local nes = r(sum)
                                drop `ww'
                            }
                        }
                        if `nes' != e(N) local rc 459
                    }
                    if !`rc' {
                        capture estimates esample: if `es', replace
                        local rc = _rc
                    }
                    if `rc' local why "sm"
                    else {
                        capture quietly stpaf `if' `in', `options' level(`level') citrans(`citrans') nodots misavedchk
                        local rc = _rc
                        if `rc' == 1 error 1
                        if !`rc' local saved 1
                        else if `rc' == 7721 local why "bs"
                        else local why "er"
                    }
                }
            }
          }
          // reason of the last failed candidate
          if !`saved' {
            if "`why'" == "ld" local ++nld
            else if "`why'" == "cm" local ++ncm
            else if "`why'" == "sm" local ++nsm
            else if "`why'" == "bs" local ++nbs
            else local ++ner
          }
        }
        local rc 0
        if `saved' local ++nsaved
        else {
            capture quietly `stcmd'
            local rc = _rc
            if !`rc' {
                capture noisily quietly stpaf `if' `in', `options' level(`level') citrans(`citrans') nodots
                local rc = _rc
            }
        }
        if !`rc' {
            if `j' == 1 matrix `BSUM' = e(b)
            else {
                capture matrix `BSUM' = `BSUM' + e(b)
                local rc = _rc
            }
        }
        if "`nodots'" == "" _dots `j' `rc'
        if `rc' {
            if "`nodots'" == "" di
            di as err "stpaf failed in imputation `m' (rc = `rc')"
            exit `rc'
        }
        matrix `TABm' = r(table)
        matrix `ALL' = nullmat(`ALL') \ `TABm'
        if `j' == 1 matrix `VPS' = r(V_paf)
        else matrix `VPS' = `VPS' + r(V_paf)
        if `j' == 1 {
            local vce0 "`r(vce)'"
            local jkm0 "`r(jk_method)'"
            local cfac "`r(cfactual)'"
            local clv "`r(clustvar)'"
            local G0 = r(N_clust)
            local Gmin = r(N_clust)
            local k = rowsof(`TABm')
            local nsc = `TABm'[`k', 1]
            forvalues s = 1/`nsc' {
                local lab`s' "`r(scen`s')'"
            }
            local lvar "`r(levelvar)'"
            local lref = r(reflevel)
            local lvnr "`r(levels)'"
            local K = r(K)
            if `K' >= . local K 0
            if `K' > 0 {
                local flist "`r(factors)'"
                matrix `WSH' = r(W_shap)
                matrix `SUBM' = r(subsets)
            }
        }
        else if "`r(vce)'" != "`vce0'" {
            if "`nodots'" == "" di
            di as err "the variance method differs between imputations (`vce0' in imputation" ///
                " `mfirst', `r(vce)' in imputation `m');"
            di as err "the default cluster jackknife probably failed in some imputations;" ///
                " specify {bf:vce(if)} or {bf:vce(jackknife)}"
            exit 498
        }
        local Gmin = min(`Gmin', r(N_clust))
        if r(jk_nrefit) < . local nref = `nref' + r(jk_nrefit)
        local nfail = `nfail' + r(N_fail)
    }
    if "`nodots'" == "" di
    restore
    _estimates unhold `H'

    // the stcox fits used must reproduce mi estimate (same data, imputations and model)
    matrix `BSUM' = `BSUM'/`M'
    if colsof(`BSUM') != colsof(`BMI') {
        di as err "the {bf:stcox} fits do not match {bf:e(b_mi)} of {bf:mi estimate}"
        exit 498
    }
    if mreldif(`BSUM', `BMI') > 1e-6 {
        di as err "the {bf:stcox} fits do not reproduce {bf:e(b_mi)} (max. rel. difference " ///
            %9.2g mreldif(`BSUM', `BMI') ");"
        di as err "the data or the imputations have changed since {bf:mi estimate}"
        exit 498
    }

    // Rubin's rules (Shapley rows on the PAF scale)
    local nuc = cond("`vce0'" == "jackknife", `Gmin' - 1, .)
    mata: _stpaf_micomb("`ALL'", `M', `nuc', `level', `logci', `K', "`VPS'", "`OUT'", "`PF'", "`SE'", "`VP'")
    // level PAFs and Shapley PAFs: largest |sum of the combined rows - combined joint PAF| over
    // the horizons; exact (up to rounding) with citrans(none), approximate with citrans(log)
    tempname SDIF SSUM
    scalar `SDIF' = .
    local ntm = rowsof(`OUT')/`nsc'
    local s1 0
    if "`lvar'" != "" {
        local s1 2
        local s2 `nsc'
    }
    if `K' > 0 {
        local s1 2
        local s2 = `K' + 1
    }
    if `s1' {
        scalar `SDIF' = 0
        forvalues h = 1/`ntm' {
            scalar `SSUM' = 0
            forvalues s = `s1'/`s2' {
                scalar `SSUM' = `SSUM' + `OUT'[(`s'-1)*`ntm' + `h', 3]
            }
            scalar `SDIF' = max(`SDIF', abs(`SSUM' - `OUT'[`h', 3]))
        }
    }

    // ---------------------------------------------------------------- display
    di
    di as txt "Population attributable fraction after {bf:mi estimate: stcox}"
    di as txt "Estimand: 1 - E[cf. events] / E[observed events]" ///
        _col(55) "Imputations" _col(70) "= " as res %9.0fc `M'
    di as txt "Counterfactual: " as res abbrev("`cfac'", 34) ///
        as txt _col(55) "Clusters" _col(70) "= " as res %9.0fc `G0'
    di as txt _col(55) "Failures (avg)" _col(70) "= " as res %9.1fc `nfail'/`M'
    if `nsaved' == `M' {
        di as txt "stcox fits: saved by {bf:mi estimate} (" as res `"`mif'"' as txt ")"
    }
    else if `nsaved' > 0 {
        di as txt "stcox fits: " as res `nsaved' as txt " saved (" as res `"`mif'"' as txt "), " ///
            as res `M' - `nsaved' as txt " refitted"
    }
    else di as txt "stcox fits: refitted in each imputation"
    if "`vce0'" == "jackknife" {
        di as txt "Within-imputation std. err.: cluster jackknife (" as res "`jkm0'" as txt ")" _c
    }
    else di as txt "Within-imputation std. err.: influence function" _c
    di as txt ", clusters: " as res cond("`clv'" == "", "record", "`clv'")
    if `nuc' < . {
        di as txt "Rubin's rules on the " cond(`logci', "log(1-PAF)", "PAF") " scale;" ///
            " df: Barnard-Rubin small sample, complete df = " `nuc'
    }
    else {
        di as txt "Rubin's rules on the " cond(`logci', "log(1-PAF)", "PAF") " scale; df: large sample"
    }
    di
    di as txt "{hline 13}{c TT}{hline 70}"
    di as txt "{ralign 12:Scenario} {c |}" _col(17) "Horizon" _col(28) "PAF" _col(37) "Std.Err." ///
        _col(48) "[`level'% Conf. Interval]" _col(76) "df" _col(81) "FMI"
    di as txt "{hline 13}{c +}{hline 70}"
    local nr = rowsof(`OUT')
    forvalues i = 1/`nr' {
        local sj = `OUT'[`i', 1]
        local h  = `OUT'[`i', 2]
        if `K' > 0 {
            if `sj' > `K' + 1 & "`subsets'" == "" continue
            if mod(`i' - 1, `ntm') == 0 {
                if `sj' == 2 di as txt "{ralign 12:Shapley} {c |}"
                else if `sj' == `K' + 2 di as txt "{ralign 12:Subsets} {c |}"
            }
        }
        local dl : subinstr local lab`sj' "Shapley: " ""
        local dfi = `OUT'[`i', 11]
        local dfs = cond(`dfi' < ., string(`dfi', "%8.1f"), "large")
        di as txt "{ralign 12:" abbrev("`dl'", 12) "} {c |}" as res _col(17) %7.4g `h' ///
            _col(26) %8.4f `OUT'[`i', 3] _col(36) %8.4f `OUT'[`i', 4] ///
            _col(48) "[" %7.4f `OUT'[`i', 5] "," %7.4f `OUT'[`i', 6] "]" ///
            _col(70) %8s "`dfs'" _col(79) %5.3f `OUT'[`i', 13]
    }
    di as txt "{hline 13}{c BT}{hline 70}"
    if "`lvar'" != "" {
        di as txt "Level PAFs ({it:k}.`lvar'): records at level {it:k} moved to `lvar' = `lref', other records"
        if `logci' {
            di as txt "unchanged.  Combined on the log(1-PAF) scale, the level PAFs add up to the joint PAF"
            di as txt "only approximately (largest difference " as res %9.2g scalar(`SDIF') as txt ");" ///
                " with {bf:citrans(none)} they add up exactly."
        }
        else di as txt "unchanged; they add up to the joint PAF (`lvar'=`lref')."
    }
    if `K' > 0 {
        di as txt "Shapley PAFs: the extra PAF from also changing the factor, averaged over all `K'! orders"
        di as txt "of changing the factors (Eide and Gefeller 1995); combined on the PAF scale."
        if `logci' {
            di as txt "The joint row is combined on the log(1-PAF) scale, so the Shapley PAFs add up to it"
            di as txt "only approximately (largest difference " as res %9.2g scalar(`SDIF') as txt ");" ///
                " with {bf:citrans(none)} they add up exactly."
        }
        else di as txt "They add up to the joint PAF."
        if "`subsets'" == "" {
            di as txt "The " 2^`K'-2 " other subset scenarios are in r(table) (option {bf:subsets} displays them)."
        }
    }
    local ciwho = cond(`K' > 0, "Joint and subset rows: CI", "CI")
    if `logci' di as txt "`ciwho' is computed on the log(1-PAF) scale and back-transformed, with t(df) quantiles."
    else di as txt "CI uses t(df) quantiles."
    if "`vce0'" == "if" & `G0' < 100 & !`vceopt' {
        di as txt "Note: the default cluster jackknife (fewer than 100 clusters) failed in the imputations;"
        di as txt "      within-imputation standard errors are from the influence function."
    }
    if `esv' {
        di as txt "Note: the estimation sample of {bf:mi estimate} varies across imputations;" ///
            " the smallest number of clusters gives the complete df."
    }
    if `nref' > 0 {
        di as txt "Note: `nref' jackknife replicate(s) refitted by {bf:stcox} after failed Newton fits (all imputations)."
    }
    // why saved fits were not used
    if `fmiss' & "`mirefit'" == "" {
        di as txt `"Note: file {bf:`mif'} of {bf:mi estimate}, {bf:saving()} not found (working directory changed?);"'
        di as txt "      {bf:stcox} was refitted in each imputation; option {bf:mifile()} gives the file's location."
    }
    if `nld' > 0 {
        di as txt `"Note: the fits of `nld' imputation(s) could not be loaded from {bf:`mif'}; {bf:stcox} was refitted there."'
    }
    if `ncm' > 0 {
        di as txt `"Note: the fits in {bf:`mif'} are not of this {bf:stcox} command (`ncm' imputation(s)); {bf:stcox} was refitted there."'
        di as txt `"      saved: `cmdx'"'
        di as txt `"      mi estimate: `stcmd'"'
    }
    if `nsm' > 0 {
        di as txt "Note: the estimation sample of the saved fit could not be reproduced in `nsm' imputation(s);" ///
            " {bf:stcox} was refitted there."
    }
    if `nbs' > 0 {
        di as txt "Note: in `nbs' imputation(s) the saved coefficients do not maximize the partial likelihood of the"
        di as txt "      data in memory (data changed since {bf:mi estimate}, or {bf:exactm}/{bf:exactp} ties with tied"
        di as txt "      failure times, which cannot be checked); {bf:stcox} was refitted there."
    }
    if `ner' > 0 {
        di as txt "Note: {bf:stpaf} failed with the saved fit in `ner' imputation(s); {bf:stcox} was refitted there."
    }

    // --------------------------------------------------------------- returns
    local last = `nr'
    return scalar paf     = `OUT'[`last', 3]
    return scalar se      = `OUT'[`last', 4]
    return scalar lb      = `OUT'[`last', 5]
    return scalar ub      = `OUT'[`last', 6]
    return scalar E       = `OUT'[`last', 7]
    return scalar E_cf    = `OUT'[`last', 8]
    return scalar se_logR = `OUT'[`last', 9]
    return scalar df      = `OUT'[`last', 11]
    return scalar riv     = `OUT'[`last', 12]
    return scalar fmi     = `OUT'[`last', 13]
    return scalar M       = `M'
    return scalar N_clust = `G0'
    return scalar N_fail  = `nfail'/`M'
    return scalar level   = `level'
    return scalar mi_nsaved = `nsaved'
    if `nuc' < . return scalar df_c = `nuc'
    if "`vce0'" == "jackknife" & "`jkm0'" != "refit" return scalar jk_nrefit = `nref'
    local vpn ""
    forvalues i = 1/`nr' {
        local vpn `vpn' s`=`OUT'[`i', 1]'h`=mod(`i' - 1, `ntm') + 1'
    }
    matrix rownames `VP' = `vpn'
    matrix colnames `VP' = `vpn'
    return matrix V_paf = `VP'
    matrix colnames `OUT' = scenario horizon paf se lb ub E E_cf se_logR events df riv fmi
    return matrix table = `OUT'
    return matrix paf_mi = `PF'
    return matrix se_mi = `SE'
    return local cfactual "`cfac'"
    return local citrans "`citrans'"
    return local vce "`vce0'"
    return local clustvar "`clv'"
    forvalues s = 1/`nsc' {
        return local scen`s' "`lab`s''"
    }
    if "`lvar'" != "" {
        return local levelvar "`lvar'"
        return scalar reflevel = `lref'
        return local levels "`lvnr'"
        return scalar sumdiff = `SDIF'
    }
    if `K' > 0 {
        return scalar K = `K'
        return local factors "`flist'"
        return matrix W_shap = `WSH'
        return matrix subsets = `SUBM'
        return scalar sumdiff = `SDIF'
    }
    if "`vce0'" == "jackknife" return local jk_method "`jkm0'"
    return local mi_fits = cond(`nsaved' == `M', "saved", cond(`nsaved' == 0, "refit", "mixed"))
    if `havef' return local mi_file `"`mif'"'
    return local m_mi "`mlist'"
    return local mi "mi"
    return local cmd "stpaf"
end

// estimation sample of a saved stcox fit, rebuilt from the stcox command line:
// _st == 1, the command's if/in, no missing values in the covariates (base variables of
// factor-variable terms), in strata() and in the vce(cluster) variable.
// Syntax: _stpaf_misamp markervar stcox <command>; markervar is set to 0/1
program define _stpaf_misamp
    version 15.0
    gettoken es 0 : 0
    gettoken cmd 0 : 0
    syntax [varlist(fv default=none)] [if] [in] [, STrata(varlist) VCE(string) *]
    marksample tm, novarlist
    quietly replace `es' = (`tm' & _st == 1)
    if "`varlist'" != "" {
        fvrevar `varlist', list
        markout `es' `r(varlist)'
    }
    if "`strata'" != "" markout `es' `strata', strok
    gettoken vt vv : vce
    local vv = strtrim("`vv'")
    if inlist("`vt'", "cl", "clu", "clus", "clust", "cluste", "cluster") & "`vv'" != "" {
        markout `es' `vv', strok
    }
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
    real colvector o1
    real colvector k1
    real colvector o0
    real colvector k0
    real matrix    in1
    real matrix    in0
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
    real scalar n, p, G, S, Ktot, nE, i, s, k, r, ns, Ks, gb, j
    real colvector eta, a, ag, ev, es, et, ew, oE, kid, newk, ks, ord, rs, o1, o0, k1, k0
    real colvector tk, dk, Kst, gbase, S0k, T0k, hk, lamk, Lc, Hc, DHc, DKc, rm0, rm1, gbr
    real colvector cp1, cp0, cr, fi, oc
    real matrix I, info, in1, in0, cinfo, RS, Eb, Ak, Gc, Dc, Ur, Ucl

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

    // ---- positions of the records on the event grid of their stratum (number of event
    //      times <= t0 and <= t1); a record is at risk at event k when rm0 < k <= rm1
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
    }
    D.S = S
    D.Ktot = Ktot
    D.Kst = Kst
    D.gbase = gbase
    gbr = gbase[strat, 1]
    // panels for the risk-set sums (0.14.0): records at risk at one or more event times, by the
    // last event time at risk (rm1) and by the last event time before entry (rm0)
    _stpaf_rpanel(selectindex(rm1 :> rm0), gbr + rm1, o1, in1, k1)
    _stpaf_rpanel(selectindex((rm1 :> rm0) :& (rm0 :> 0)), gbr + rm0, o0, in0, k0)
    D.o1 = o1
    D.in1 = in1
    D.k1 = k1
    D.o0 = o0
    D.in0 = in0
    D.k0 = k0

    // ---- risk-set sums at every event time (vectorized, 0.14.0): S0, S1 over all records,
    //      T0, T1 over the records of the if subset (g); the counterfactual sums of _stpaf_scen
    //      use the same arithmetic, so a scenario that changes nothing gives E0 = E exactly
    ag = a :* g
    RS = _stpaf_rsum(D, (a, a :* X, ag, ag :* X))
    S0k = RS[., 1]
    Eb = RS[., 2..(p+1)] :/ S0k
    T0k = RS[., p+2]
    Ak = RS[., (p+3)..(2*p+2)] :/ S0k
    RS = J(0, 0, .)

    // ---- cumulative arrays, one block of length K_s+1 per stratum;
    //      element m (0..K_s) of stratum s lives at index gbase[s] + s + m
    lamk = dk :/ S0k
    hk = T0k :/ S0k
    Lc = _stpaf_cum0(D, lamk)
    Hc = _stpaf_cum0(D, hk :* lamk)
    DHc = _stpaf_cum0(D, dk :* hk)
    DKc = _stpaf_cum0(D, dk)
    Gc = J(Ktot + S, p, 0)
    Dc = J(Ktot + S, p, 0)
    for (j = 1; j <= p; j++) {
        Gc[., j] = _stpaf_cum0(D, lamk :* Eb[., j])
        Dc[., j] = _stpaf_cum0(D, dk :* (Ak[., j] - hk :* Eb[., j]))
    }

    // ---- information and cluster-level partial-likelihood score (Breslow):
    //      I = X' diag(a c) X - sum_k dk Eb Eb', c_r = Lc summed over the events where r is at risk
    cp1 = gbr + strat + rm1
    cp0 = gbr + strat + rm0
    cr = Lc[cp1, 1] - Lc[cp0, 1]
    I = quadcross(X, a :* cr, X) - quadcross(Eb, dk, Eb)
    Ur = -a :* (X :* cr - (Gc[cp1, .] - Gc[cp0, .]))
    fi = selectindex((d :!= 0) :& (w :> 0))
    if (length(fi)) Ur[fi, .] = Ur[fi, .] + w[fi, 1] :* (X[fi, .] - Eb[gbr[fi, 1] + rm1[fi, 1], .])
    oc = order(clus, 1)
    cinfo = panelsetup(clus[oc, 1], 1)
    Ucl = J(G, p, 0)
    Ucl[clus[oc[cinfo[., 1], 1], 1], .] = panelsum(Ur[oc, .], cinfo)
    Ur = J(0, 0, .)

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

// ------------------------------------------- vectorized risk-set sums (0.14.0)
// records idx sorted by key (event index of the stratum's grid); o: sorted records, info:
// panels of equal key, k: the key of each panel
void _stpaf_rpanel(real colvector idx, real colvector key, real colvector o, real matrix info,
    real colvector k)
{
    real colvector kk, ord, ks

    if (length(idx) == 0) {
        o = J(0, 1, .)
        info = J(0, 2, .)
        k = J(0, 1, .)
        return
    }
    kk = key[idx, 1]
    ord = order(kk, 1)
    o = idx[ord, 1]
    ks = kk[ord, 1]
    info = panelsetup(ks, 1)
    k = ks[info[., 1], 1]
}

// risk-set sums at every event time (records with rm0 < k <= rm1) of the columns of Y
// (records x columns): bucket sums by the last event at risk minus bucket sums by the last
// event before entry, cumulated backwards within strata (as _stpaf_jkrisk)
real matrix _stpaf_rsum(struct stpaf_dat scalar D, real matrix Y)
{
    real matrix B, out
    real colvector v, r
    real scalar s, K, lo, hi, j

    B = _stpaf_bsum(Y[D.o1, .], D.in1, D.k1, D.Ktot)
    if (rows(D.in0)) B = B - _stpaf_bsum(Y[D.o0, .], D.in0, D.k0, D.Ktot)
    out = B
    for (s = 1; s <= D.S; s++) {
        K = D.Kst[s]
        if (K == 0) continue
        lo = D.gbase[s] + 1
        hi = D.gbase[s] + K
        for (j = 1; j <= cols(B); j++) {
            v = B[|lo, j \ hi, j|]
            r = quadrunningsum(v[(K::1), 1])
            out[|lo, j \ hi, j|] = r[(K::1), 1]
        }
    }
    return(out)
}

// forward cumulation of v (one value per event time) with a leading zero per stratum
// (layout of D.Lc)
real colvector _stpaf_cum0(struct stpaf_dat scalar D, real colvector v)
{
    real colvector out
    real scalar s, K, pb

    out = J(D.Ktot + D.S, 1, 0)
    for (s = 1; s <= D.S; s++) {
        K = D.Kst[s]
        if (K == 0) continue
        pb = D.gbase[s] + s
        out[|pb + 1 \ pb + K|] = quadrunningsum(v[|D.gbase[s] + 1 \ D.gbase[s] + K|])
    }
    return(out)
}

// ------------------------------------------- counterfactual columns (0.13.0)
// Plan of the column versions. fvl: the non-omitted e(b) terms (columns of X); svs: the
// (unabbreviated) variable of each change i = 1..nspec; Stata locals scl1..scl<nsc>: the changes
// of each scenario. Column j depends on the changes C_j of the variables in its term (the name
// after the last "." of each "#" part); in scenario s it is the version (j, S_s & C_j), or the
// observed column when that set is empty. Versions are numbered in order of first appearance.
// Rebuild passes: a pass with the changes of scenario s gives the version (j, S_s & C_j) of
// every column; passes are chosen greedily (the scenario giving most missing versions, ties to
// the first scenario) until all versions are covered. Results: externals __stpaf_XV (n x nu,
// filled by _stpaf_xcput), __stpaf_XM (nsc x p version indices), __stpaf_XC (1 x nu column of
// each version); locals pfx+"nu", pfx+"nb" and for pass b: pfx+"s"+b (scenario), pfx+"c"+b
// (columns to store), pfx+"v"+b (their version numbers), pfx+"k"+b (columns that must equal
// the observed ones). Checked in Python (xcmem_check.py).
void _stpaf_xcplan(string scalar fvl, string scalar svs, real scalar nsc, real scalar n,
    string scalar pfx)
{
    transmorphic scalar A
    pointer(real matrix) scalar pV, pM, pC
    string rowvector terms, sv, parts
    string scalar key, v
    real matrix CJ, SCM, KI
    real rowvector idx, T, cov, used, nz, newc, newv, chk, vcol
    real scalar p, ns, s, j, i, k, nu, nb, best, bs, gsum

    terms = tokens(fvl)
    p = cols(terms)
    sv = tokens(svs)
    ns = cols(sv)
    CJ = J(ns, p, 0)
    for (j = 1; j <= p; j++) {
        parts = tokens(subinstr(terms[j], "#", " ", .))
        for (k = 1; k <= cols(parts); k++) {
            v = regexr(parts[k], "^.*[.]", "")
            for (i = 1; i <= ns; i++) {
                if (sv[i] == v) CJ[i, j] = 1
            }
        }
    }
    SCM = J(nsc, ns, 0)
    for (s = 1; s <= nsc; s++) {
        idx = strtoreal(tokens(st_local("scl" + strofreal(s))))
        if (cols(idx)) SCM[s, idx] = J(1, cols(idx), 1)
    }
    // unique (column, set of changes) versions
    A = asarray_create()
    KI = J(nsc, p, 0)
    vcol = J(1, 0, .)
    nu = 0
    for (s = 1; s <= nsc; s++) {
        for (j = 1; j <= p; j++) {
            T = SCM[s, .] :* (CJ[., j]')
            if (!any(T)) continue
            key = strofreal(j) + ":" + invtokens(strofreal(selectindex(T)))
            if (asarray_contains(A, key)) KI[s, j] = asarray(A, key)
            else {
                nu++
                asarray(A, key, nu)
                KI[s, j] = nu
                vcol = vcol, j
            }
        }
    }
    // rebuild passes
    cov = J(1, nu, 0)
    used = J(1, nsc, 0)
    nb = 0
    while (nu > 0) {
        best = 0
        bs = 0
        for (s = 1; s <= nsc; s++) {
            if (used[s]) continue
            nz = selectindex(KI[s, .] :> 0)
            if (length(nz) == 0) continue
            gsum = sum(cov[KI[s, nz]] :== 0)
            if (gsum > best) {
                best = gsum
                bs = s
            }
        }
        if (best == 0) break
        used[bs] = 1
        nb++
        nz = selectindex(KI[bs, .] :> 0)
        newc = select(nz, cov[KI[bs, nz]] :== 0)
        newv = KI[bs, newc]
        cov[newv] = J(1, length(newv), 1)
        chk = selectindex(KI[bs, .] :== 0)
        st_local(pfx + "s" + strofreal(nb), strofreal(bs))
        st_local(pfx + "c" + strofreal(nb), invtokens(strofreal(newc)))
        st_local(pfx + "v" + strofreal(nb), invtokens(strofreal(newv)))
        st_local(pfx + "k" + strofreal(nb), (length(chk) ? invtokens(strofreal(chk)) : ""))
    }
    pV = findexternal("__stpaf_XV")
    if (pV == NULL) pV = crexternal("__stpaf_XV")
    pM = findexternal("__stpaf_XM")
    if (pM == NULL) pM = crexternal("__stpaf_XM")
    pC = findexternal("__stpaf_XC")
    if (pC == NULL) pC = crexternal("__stpaf_XC")
    *pV = J(n, nu, .)
    *pM = KI
    *pC = vcol
    st_local(pfx + "nu", strofreal(nu))
    st_local(pfx + "nb", strofreal(nb))
}

// one rebuild pass: store the columns vnames (versions vids) in __stpaf_XV, and compare the
// columns cnames (positions cpos) with the observed design matrix __stpaf_X0; chkname: 0, or the
// position of the first column that differs (internal error)
void _stpaf_xcput(string scalar vnames, string scalar vids, string scalar cnames,
    string scalar cpos, string scalar sel, string scalar chkname)
{
    pointer(real matrix) scalar pV, pX
    real rowvector v, c, dd, bad
    real matrix Y

    pV = findexternal("__stpaf_XV")
    pX = findexternal("__stpaf_X0")
    if (pV == NULL | pX == NULL) _error("stpaf internal error: counterfactual columns not set up")
    if (vnames != "") {
        v = strtoreal(tokens(vids))
        (*pV)[., v] = st_data(., tokens(vnames), sel)
    }
    st_numscalar(chkname, 0)
    if (cnames != "") {
        c = strtoreal(tokens(cpos))
        Y = st_data(., tokens(cnames), sel)
        dd = colsum(Y :!= (*pX)[., c])
        if (any(dd)) {
            bad = selectindex(dd :> 0)
            st_numscalar(chkname, c[bad[1]])
        }
    }
}

// counterfactual design matrix of scenario s: X with the stored versions in their columns
real matrix _stpaf_xcget(real matrix X, real scalar s)
{
    pointer(real matrix) scalar pV, pM
    real rowvector m, c
    real matrix Xc

    pV = findexternal("__stpaf_XV")
    pM = findexternal("__stpaf_XM")
    if (pV == NULL | pM == NULL) _error("stpaf internal error: counterfactual columns not found")
    m = (*pM)[s, .]
    Xc = X
    c = selectindex(m :> 0)
    if (length(c)) Xc[., c] = (*pV)[., m[c]]
    return(Xc)
}

// PAF (one row per horizon) for the counterfactual design matrix Xc
// columns: horizon, paf, se, lb, ub, E, E0, se(log R), failures up to horizon; se, lb, ub and
// se(log R) are left missing (filled from the covariance by _stpaf_run). Outputs for the
// influence function (_stpaf_ifv): h0k (counterfactual h at the event times), H0c (cumulative
// sum of h0k dL, layout of D.Lc) and D0t (horizons x p: dE0/db up to each horizon)
real matrix _stpaf_scen(struct stpaf_dat scalar D, real matrix Xc, real rowvector tau,
    real colvector h0k, real colvector H0c, real matrix D0t)
{
    real scalar p, S, Ktot, nt, s, j, Ks, gb, pb, ti, mts, E, E0, evn
    real colvector ac, lamk, DH0c
    real rowvector d0
    real matrix RS, Ack, D0c, res

    p = D.p
    S = D.S
    Ktot = D.Ktot
    nt = cols(tau)

    ac = D.w :* D.g :* exp(Xc*D.b' :- D.cen)

    // ---- counterfactual risk-set sums at every event time (vectorized, 0.14.0; same
    //      arithmetic as the observed sums of _stpaf_prep)
    RS = _stpaf_rsum(D, (ac, ac :* Xc))
    h0k = RS[., 1] :/ D.S0k
    Ack = RS[., 2..(p+1)] :/ D.S0k
    RS = J(0, 0, .)
    lamk = D.dk :/ D.S0k
    H0c = _stpaf_cum0(D, h0k :* lamk)
    DH0c = _stpaf_cum0(D, D.dk :* h0k)
    D0c = J(Ktot + S, p, 0)
    for (j = 1; j <= p; j++) D0c[., j] = _stpaf_cum0(D, D.dk :* (Ack[., j] - h0k :* D.Eb[., j]))

    // ---- point estimate for each horizon
    res = J(nt, 9, .)
    D0t = J(nt, p, 0)
    for (ti = 1; ti <= nt; ti++) {
        E = 0
        E0 = 0
        evn = 0
        d0 = J(1, p, 0)
        for (s = 1; s <= S; s++) {
            Ks = D.Kst[s]
            gb = D.gbase[s]
            pb = gb + s
            mts = 0
            if (Ks > 0) mts = _stpaf_pos(D.tk, gb, Ks, tau[ti])
            E = E + D.DHc[pb + mts]
            E0 = E0 + DH0c[pb + mts]
            evn = evn + D.DKc[pb + mts]
            d0 = d0 + D0c[pb + mts,.]
        }
        if (E <= 0 | E0 <= 0) _error("no expected events in the selected subset/horizon")
        D0t[ti, .] = d0
        res[ti,1] = tau[ti]
        if (tau[ti] >= .) res[ti,1] = max(D.tk)
        res[ti,2] = 1 - E0/E
        res[ti,6] = E
        res[ti,7] = E0
        res[ti,9] = evn
    }
    return(res)
}

// covariance of the PAFs of all scenarios and horizons (columns (s-1)*nt + t) from the cluster
// influence functions Psi = -(1 - PAF) phi, phi = psi0/E0 - psi/E (0.13.0: replaces the record
// loop of _stpaf_scen and the clusters x columns matrix of 0.12.0). The clusters are processed
// in blocks of chunk clusters (default: about 2^24 elements of Psi per block); for each block
// the influence functions of all columns are built with vectorized expressions and
// Psi_block'Psi_block is accumulated. V = G/(G-1) Psi'Psi (fwrep: N/(N-1) Psi' diag(1/w) Psi
// over w > 0, one record per cluster). Per record r (stratum s, positions q0, q1 of its interval
// truncated at the horizon, in the cumulative arrays):
//   v  = -a_r (Hc[q1] - Hc[q0]) + g_r a_r (Lc[q1] - Lc[q0]) + [failure up to the horizon] w_r hk
//   v0 = -a_r (H0c[q1] - H0c[q0]) + ac_r (Lc[q1] - Lc[q0]) + [failure up to the horizon] w_r h0k
// psi = cluster sum of v + IFb Dt', psi0 likewise with D0t. ac_r = w g exp(eta^c - cen) with
// eta^c = X b + (versions - X) b over the changed columns. H0C, H0K: one column per scenario;
// D0T, E0V: one row per column of Psi; EV: E per horizon. nch: number of blocks.
// Checked in Python against the 0.12.0 loop (xcmem_check.py).
real matrix _stpaf_ifv(struct stpaf_dat scalar D, real matrix X, real rowvector tau,
    real scalar nsc, real matrix H0C, real matrix H0K, real matrix D0T, real colvector EV,
    real colvector E0V, real scalar fwrep, real scalar chunk, real scalar nch)
{
    pointer(real matrix) scalar pV, pM
    real matrix mt, Dt, V, PS, PSO, P1M, P0M, DLM, FM, info, infc, IFc
    real colvector ordc, wcl, rr, sr, pbr, ar, gr, wr, m0r, m1r, fr, FKr, eta, etac, ac
    real colvector mtr, q1, q0, P1, P0, dL, Fi, iF, vo, v0, h0, psi0, wc, ii
    real rowvector cs
    real scalar G, S, p, nt, m, nc, gc, s, t, g1, g2, r1, r2, col, R, Nw, cks

    G = D.G
    S = D.S
    p = D.p
    nt = cols(tau)
    m = nsc*nt
    mt = J(S, nt, 0)
    for (s = 1; s <= S; s++) {
        if (D.Kst[s] == 0) continue
        for (t = 1; t <= nt; t++) mt[s, t] = _stpaf_pos(D.tk, D.gbase[s], D.Kst[s], tau[t])
    }
    Dt = J(nt, p, 0)
    for (t = 1; t <= nt; t++) {
        for (s = 1; s <= S; s++) Dt[t, .] = Dt[t, .] + D.Dc[D.gbase[s] + s + mt[s, t], .]
    }
    pV = findexternal("__stpaf_XV")
    pM = findexternal("__stpaf_XM")
    if (pV == NULL | pM == NULL) _error("stpaf internal error: counterfactual columns not found")
    ordc = order(D.clus, 1)
    info = panelsetup(D.clus[ordc, 1], 1)
    if (rows(info) != G) _error("stpaf internal error: cluster codes are not 1..G")
    if (fwrep) {
        wcl = J(G, 1, 0)
        wcl[D.clus, 1] = D.w
    }
    cks = (chunk >= 1 ? chunk : max((1, floor(2^24/m))))
    V = J(m, m, 0)
    nch = 0
    for (g1 = 1; g1 <= G; g1 = g2 + 1) {
        g2 = min((G, g1 + cks - 1))
        nch++
        r1 = info[g1, 1]
        r2 = info[g2, 2]
        rr = ordc[|r1 \ r2|]
        infc = info[|g1, 1 \ g2, 2|] :- (r1 - 1)
        gc = g2 - g1 + 1
        nc = rows(rr)
        sr = D.strat[rr, 1]
        pbr = D.gbase[sr, 1] + sr
        ar = D.a[rr, 1]
        gr = D.g[rr, 1]
        wr = D.w[rr, 1]
        m0r = D.rm0[rr, 1]
        m1r = D.rm1[rr, 1]
        fr = (D.d[rr, 1] :!= 0) :& (wr :> 0)
        FKr = D.gbase[sr, 1] + m1r
        eta = D.xb[rr, 1]
        IFc = D.IFb[|g1, 1 \ g2, p|]
        // observed side and record positions, per horizon
        P1M = J(nc, nt, .)
        P0M = J(nc, nt, .)
        DLM = J(nc, nt, .)
        FM = J(nc, nt, 0)
        PSO = J(gc, nt, .)
        for (t = 1; t <= nt; t++) {
            mtr = mt[sr, t]
            q1 = rowmin((m1r, mtr))
            q0 = rowmin((m0r, mtr))
            P1 = pbr + q1
            P0 = pbr + q0
            dL = D.Lc[P1, 1] - D.Lc[P0, 1]
            Fi = fr :& (m1r :<= mtr)
            iF = selectindex(Fi)
            vo = -ar :* (D.Hc[P1, 1] - D.Hc[P0, 1]) + gr :* ar :* dL
            if (length(iF)) vo[iF, 1] = vo[iF, 1] + wr[iF, 1] :* D.hk[FKr[iF, 1], 1]
            PSO[., t] = panelsum(vo, infc) + IFc*Dt[t, .]'
            P1M[., t] = P1
            P0M[., t] = P0
            DLM[., t] = dL
            FM[., t] = Fi
        }
        // counterfactual side, every scenario and horizon
        PS = J(gc, m, .)
        for (s = 1; s <= nsc; s++) {
            cs = selectindex((*pM)[s, .] :> 0)
            if (length(cs)) etac = eta + ((*pV)[rr, (*pM)[s, cs]] - X[rr, cs])*D.b[1, cs]'
            else etac = eta
            ac = wr :* gr :* exp(etac :- D.cen)
            h0 = H0C[., s]
            for (t = 1; t <= nt; t++) {
                col = (s - 1)*nt + t
                v0 = -ar :* (h0[P1M[., t], 1] - h0[P0M[., t], 1]) + ac :* DLM[., t]
                iF = selectindex(FM[., t])
                if (length(iF)) v0[iF, 1] = v0[iF, 1] + wr[iF, 1] :* H0K[FKr[iF, 1], s]
                psi0 = panelsum(v0, infc) + IFc*D0T[col, .]'
                R = E0V[col]/EV[t]
                PS[., col] = -R :* (psi0 :/ E0V[col] - PSO[., t] :/ EV[t])
            }
        }
        if (fwrep) {
            wc = wcl[|g1 \ g2|]
            ii = selectindex(wc :> 0)
            if (length(ii)) V = V + quadcross(PS[ii, .], 1 :/ wc[ii, 1], PS[ii, .])
        }
        else V = V + quadcross(PS, PS)
    }
    if (fwrep) {
        Nw = sum(D.w)
        V = Nw/(Nw - 1)*V
    }
    else V = G/(G - 1)*V
    return(V)
}

// driver: prepare once, evaluate every scenario, post results to Stata
// (the structure never leaves Mata function scope)
// tiemode 1 (efron) with tied failure times and unit weights: Efron coefficient term
// (iftname = 1); vbname: cluster-robust VCE of b, G/(G-1) IFb'IFb (not with fwrep);
// vpname: covariance of the PAFs of all scenarios and horizons (rows of TAB), from the cluster
// influence functions (_stpaf_ifv, blocks of chunk clusters; nchname: number of blocks); its
// diagonal gives se, and se(log R) = se/(1 - PAF)
void _stpaf_run(real matrix X, real colvector t0, real colvector t1, real colvector d,
    real colvector g, real colvector w, real colvector clus, real colvector strat,
    real rowvector b, real scalar nsc, real rowvector tau, real scalar level,
    real scalar logci, real scalar fwrep, string scalar tabname, string scalar scname,
    string scalar stname, string scalar chkname, string scalar xbvar, string scalar sel,
    real scalar tiemode, string scalar vbname, string scalar iftname, string scalar vpname,
    real scalar chunk, string scalar nchname)
{
    struct stpaf_dat scalar D
    real matrix res, out, H0C, H0K, D0T, D0t, VP
    real colvector h0k, H0c, rr
    real scalar j, efr, nt, m, c, z, se, R, selr, paf, nch

    D = _stpaf_prep(X, t0, t1, d, g, w, clus, strat, b)
    efr = 0
    if (tiemode == 1 & all(w :== 1) & max(D.dk) > 1) {
        _stpaf_efronif(X, D)
        efr = 1
    }
    st_numscalar(iftname, efr)
    if (!fwrep) st_matrix(vbname, D.G/(D.G - 1)*quadcross(D.IFb, D.IFb))
    st_numscalar(chkname, _stpaf_chk(D, xbvar, sel))
    _stpaf_post(D, scname, stname)
    nt = cols(tau)
    m = nsc*nt
    out = J(m, 10, .)
    H0C = J(D.Ktot + D.S, nsc, .)
    H0K = J(D.Ktot, nsc, .)
    D0T = J(m, D.p, .)
    for (j = 1; j <= nsc; j++) {
        res = _stpaf_scen(D, _stpaf_xcget(X, j), tau, h0k, H0c, D0t)
        rr = ((j - 1)*nt + 1)::(j*nt)
        out[rr, .] = (J(nt, 1, j), res)
        H0C[., j] = H0c
        H0K[., j] = h0k
        D0T[rr, .] = D0t
    }
    nch = 0
    VP = _stpaf_ifv(D, X, tau, nsc, H0C, H0K, D0T, out[1::nt, 7], out[., 8], fwrep, chunk, nch)
    z = invnormal(1 - (1 - level/100)/2)
    for (c = 1; c <= m; c++) {
        se = sqrt(VP[c, c])
        R = out[c, 8]/out[c, 7]
        paf = out[c, 3]
        selr = se/R
        out[c, 4] = se
        out[c, 9] = selr
        if (logci) {
            out[c, 5] = 1 - R*exp(z*selr)
            out[c, 6] = 1 - R*exp(-z*selr)
        }
        else {
            out[c, 5] = paf - z*se
            out[c, 6] = paf + z*se
        }
    }
    st_numscalar(nchname, nch)
    st_matrix(tabname, out)
    st_matrix(vpname, VP)
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
    real matrix PF, res, Xc, D0t
    real colvector ii, h0k, H0c
    real scalar nt, i, gg, j

    nt = cols(tau)
    PF = J(rows(cl), nsc*nt, .)
    for (i = 1; i <= rows(cl); i++) {
        gg = cl[i]
        ii = selectindex(clus :!= gg)
        D = _stpaf_prep(X[ii,.], t0[ii], t1[ii], d[ii], g[ii], w[ii],
            _stpaf_recode(clus[ii]), _stpaf_recode(strat[ii]), BJ[i,.])
        for (j = 1; j <= nsc; j++) {
            Xc = _stpaf_xcget(X, j)
            res = _stpaf_scen(D, Xc[ii, .], tau, h0k, H0c, D0t)
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
    real scalar logci, string scalar tabname, string scalar vpname)
{
    real matrix PF

    PF = _stpaf_jkrefit(X, t0, t1, d, g, w, clus, strat, BJ, (1::rows(BJ)), nsc, tau,
        level, logci)
    _stpaf_jkvar(PF, level, logci, tabname, vpname)
}

// fast jackknife with failed Newton fits: the replicate PAFs of the other clusters are in
// the external __stpaf_PF; those of the clusters cl come from the stcox refits BJ
void _stpaf_jackfill(real matrix X, real colvector t0, real colvector t1, real colvector d,
    real colvector g, real colvector w, real colvector clus, real colvector strat,
    real matrix BJ, real colvector cl, real scalar nsc, real rowvector tau, real scalar level,
    real scalar logci, string scalar tabname, string scalar vpname)
{
    pointer() scalar pPF
    real matrix PF

    pPF = findexternal("__stpaf_PF")
    if (pPF == NULL) _error("stpaf internal error: replicate PAFs of the Newton jackknife not found")
    PF = *pPF
    if (rows(BJ) != rows(cl)) _error("stpaf internal error: refitted clusters do not match")
    PF[cl, .] = _stpaf_jkrefit(X, t0, t1, d, g, w, clus, strat, BJ, cl, nsc, tau, level, logci)
    _stpaf_jkvar(PF, level, logci, tabname, vpname)
}

// jackknife variance from the replicate PAFs (G x (scenarios*horizons), in the
// row order of TAB): (G-1)/G sum (theta_-g - mean)^2 of log(1-PAF) and of PAF;
// confidence limits with the t quantile on G-1 degrees of freedom (as Stata's jackknife);
// vpname: jackknife covariance of the PAFs of all rows, (G-1)/G sum (theta - mean)(theta - mean)'
void _stpaf_jkvar(real matrix PF, real scalar level, real scalar logci, string scalar tabname,
    string scalar vpname)
{
    real matrix T
    real colvector lr, pf
    real rowvector pm
    real scalar G, c, z, R, selr, sep, paf

    G = rows(PF)
    if (hasmissing(PF)) _error("jackknife replicate PAF could not be computed")
    if (max(PF) >= 1) _error("jackknife replicate PAF >= 1; log(1-PAF) undefined")
    T = st_matrix(tabname)
    z = invttail(G - 1, (1 - level/100)/2)
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
    pm = mean(PF)
    st_matrix(tabname, T)
    st_matrix(vpname, (G-1)/G*quadcrossdev(PF, pm, PF, pm))
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
    real scalar    efron
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
    real colvector kfr
    real colvector fpan
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
    real colvector act, idx, o, k, fs
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
    // Efron: event index and panel number of every failure row (rows sorted as of)
    JK.efron = 0
    JK.kfr = JK.gb[JK.of, 1] + D.rm1[JK.of, 1]
    fs = J(rows(JK.of), 1, 0)
    if (rows(JK.inf)) fs[JK.inf[., 1], 1] = J(rows(JK.inf), 1, 1)
    JK.fpan = runningsum(fs)
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

// Efron-exact coefficient part of the influence function (0.9.0), for unit weights.
// Replaces D.B, D.IFb and D.scorechk (computed by _stpaf_prep from the Breslow score and
// information) by the Efron ones. Per record r, with the notation of _stpaf_jkllE
// (den_l, E_l = (S1 - (l/d) SF1)/den_l at an event with d failures) and
//   lamE = sum_l 1/den_l, mu = sum_l (l/d)/den_l, GE = sum_l E_l/den_l,
//   GmuE = sum_l (l/d) E_l/den_l, Ebar = mean_l E_l,
// the centred (Therneau) score residual is
//   U_r = -a_r (x_r cE_r - cGE_r) + d_r [(x_r - Ebar) + a_r (x_r mu - GmuE)] (at r's failure)
// where cE_r, cGE_r cumulate lamE, GE over the events at which r is at risk. The residuals sum
// to the Efron score, and without tied failure times they are the Breslow ones.
// Information: X' diag(a cE) X - sum_F a mu x x' - sum_l E_l E_l'.
// Checked in Python (efron_if_check.py): residuals = record-by-record loop (2e-15), sum = Efron
// score, information = numerical derivative of the score, I^-1 U_g = leave-one-cluster-out
// change in b up to O(1/G).
void _stpaf_efronif(real matrix X, struct stpaf_dat scalar D)
{
    struct stpaf_jk scalar JK
    real colvector a, a1, a0, af, dk, S0, SF0, cx, st, l, dkr, cl, den, ide, lamE, mu, cE, ce
    real colvector muf, o
    real matrix S1, SF1, V0, Er, GE, GmuE, Ebar, cGE, Ur, I, info, Ucl
    real scalar j

    JK = _stpaf_jkset(X, D)
    // a = exp(X b - cen) (w = 1); failure rows: JK.of (all with w = 1)
    a = D.a
    af = a[JK.of, 1]
    dk = _stpaf_bsum(J(rows(JK.of), 1, 1), JK.inf, JK.kf, JK.Ktot)
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
    SF0 = _stpaf_bsum(af, JK.inf, JK.kf, JK.Ktot)
    SF1 = _stpaf_bsum(af :* JK.Xf, JK.inf, JK.kf, JK.Ktot)
    // rank l (0..d-1) of each failure row within its event time
    cx = runningsum(J(rows(JK.of), 1, 1)) :- 1
    st = cx[JK.inf[., 1], 1]
    l = cx - st[JK.fpan, 1]
    dkr = dk[JK.kfr, 1]
    cl = l :/ dkr
    den = S0[JK.kfr, 1] - cl :* SF0[JK.kfr, 1]
    Er = (S1[JK.kfr, .] - cl :* SF1[JK.kfr, .]) :/ den
    ide = 1 :/ den
    lamE = _stpaf_bsum(ide, JK.inf, JK.kf, JK.Ktot)
    mu = _stpaf_bsum(ide :* cl, JK.inf, JK.kf, JK.Ktot)
    GE = _stpaf_bsum(Er :* ide, JK.inf, JK.kf, JK.Ktot)
    GmuE = _stpaf_bsum(Er :* (ide :* cl), JK.inf, JK.kf, JK.Ktot)
    Ebar = _stpaf_bsum(Er, JK.inf, JK.kf, JK.Ktot) :/ (dk + (dk :== 0))
    cE = _stpaf_jkcum0(JK, lamE)
    ce = cE[JK.cp1, 1] - cE[JK.cp0, 1]
    cGE = J(JK.Ktot + JK.S, JK.p, 0)
    for (j = 1; j <= JK.p; j++) cGE[., j] = _stpaf_jkcum0(JK, GE[., j])
    Ur = -a :* (X :* ce - (cGE[JK.cp1, .] - cGE[JK.cp0, .]))
    muf = af :* mu[JK.kfr, 1]
    Ur[JK.of, .] = Ur[JK.of, .] + (JK.Xf - Ebar[JK.kfr, .]) + (muf :* JK.Xf - af :* GmuE[JK.kfr, .])
    I = quadcross(X, a :* ce, X) - quadcross(JK.Xf, muf, JK.Xf) - quadcross(Er, Er)
    // cluster sums (codes 1..G, all present)
    o = order(D.clus, 1)
    info = panelsetup(D.clus[o, 1], 1)
    Ucl = J(D.G, D.p, 0)
    Ucl[D.clus[o[info[., 1], 1], 1], .] = panelsum(Ur[o, .], info)
    D.B = invsym(I)
    D.IFb = Ucl*D.B
    D.scorechk = abs(colsum(Ucl)) :/ sqrt(diagonal(I))'
}

// Efron log partial likelihood (up to a constant), score and information at b for an
// unweighted model (w = 1, or 0 for the left-out cluster), and the Breslow cumulative baseline
// hazard Lc, which the replicate PAF uses as the point estimate and the refit path do.
// At an event time with d failures (set F, SF0 = sum of a over F, SF1 = sum of a x over F) the
// l-th failure (l = 0..d-1) has denominator S0 - (l/d) SF0. With
// lamE = sum_l 1/den_l, mu = sum_l (l/d)/den_l and E_l = (S1 - (l/d) SF1)/den_l:
//   ll = sum_F eta - sum_l ln(den_l)
//   U  = sum_F x - sum_r a_r x_r cE_r + sum_F a x mu
//   I  = X' diag(a cE) X - sum_F a mu x x' - sum_l E_l E_l'
// cE_r: lamE cumulated over the events at which record r is at risk (as c for Breslow).
// With d = 1 everywhere this is the Breslow likelihood.
void _stpaf_jkllE(struct stpaf_jk scalar JK, real matrix X, real colvector w, real rowvector b,
    real scalar ll, real rowvector U, real matrix I, real colvector Lc)
{
    real colvector eta, a, a1, a0, wf, af, dk, S0, S0s, SF0, cx, st, l, dkr, cl, act, den
    real colvector ide, lamE, mu, cE, ce, ace, muf, ii, pos
    real matrix S1, SF1, V0, Er

    eta = X*b' :- JK.cen
    a = w :* exp(eta)
    wf = w[JK.of, 1]
    af = a[JK.of, 1]
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
    S0s = S0 + (1 :- pos)
    Lc = _stpaf_jkcum0(JK, dk :/ S0s)
    SF0 = _stpaf_bsum(af, JK.inf, JK.kf, JK.Ktot)
    SF1 = _stpaf_bsum(af :* JK.Xf, JK.inf, JK.kf, JK.Ktot)
    // rank l of each remaining failure among the remaining failures at its event time
    cx = runningsum(wf) - wf
    st = cx[JK.inf[., 1], 1]
    l = cx - st[JK.fpan, 1]
    dkr = dk[JK.kfr, 1]
    cl = l :/ (dkr + (dkr :== 0))
    act = (wf :> 0)
    // failures of the left-out cluster: no term (den set to 1, weight 0)
    den = act :* (S0[JK.kfr, 1] - cl :* SF0[JK.kfr, 1]) + (1 :- act)
    ide = wf :/ den
    lamE = _stpaf_bsum(ide, JK.inf, JK.kf, JK.Ktot)
    mu = _stpaf_bsum(ide :* cl, JK.inf, JK.kf, JK.Ktot)
    cE = _stpaf_jkcum0(JK, lamE)
    ce = cE[JK.cp1, 1] - cE[JK.cp0, 1]
    ace = a :* ce
    muf = af :* mu[JK.kfr, 1]
    Er = (S1[JK.kfr, .] - cl :* SF1[JK.kfr, .]) :/ den
    ii = selectindex(act)
    ll = quadsum(wf :* eta[JK.of, 1]) - quadsum(ln(den[ii, 1]))
    U = quadcolsum(wf :* JK.Xf) - quadcross(ace, X) + quadcross(muf, JK.Xf)
    I = quadcross(X, ace, X) - quadcross(JK.Xf, muf, JK.Xf) - quadcross(Er, wf, Er)
}

// weighted Breslow log partial likelihood (up to a constant), score, information, and
// the cumulative baseline hazard Lc, at b, with record weights w (JK.efron: Efron likelihood)
void _stpaf_jkll(struct stpaf_jk scalar JK, real matrix X, real colvector w, real rowvector b,
    real scalar ll, real rowvector U, real matrix I, real colvector Lc)
{
    real colvector eta, a, a1, a0, wf, dk, S0, S0s, pos, c, ac, ii
    real matrix S1, Eb, V0

    if (JK.efron) {
        _stpaf_jkllE(JK, X, w, b, ll, U, I, Lc)
        return
    }
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
// Q1, Q0: positions in Lc of each record's interval truncated at each horizon;
// counterfactual linear predictor (0.13.0): X b + DX[., versions] b[columns], DX = versions -
// observed columns (set up once by _stpaf_jackfast), M = __stpaf_XM
real rowvector _stpaf_jkpaf(struct stpaf_jk scalar JK, real matrix X, real colvector w,
    real colvector g, real rowvector b, real colvector Lc, real scalar nsc,
    real matrix Q1, real matrix Q0, real matrix DX, real matrix M)
{
    real rowvector out, cs
    real colvector eta, etac, a, ac, dL
    real scalar j, t, nt, E, E0

    nt = cols(Q1)
    out = J(1, nsc*nt, .)
    eta = X*b'
    a = w :* g :* exp(eta :- JK.cen)
    for (j = 1; j <= nsc; j++) {
        cs = selectindex(M[j, .] :> 0)
        if (length(cs)) etac = eta + DX[., M[j, cs]]*b[1, cs]'
        else etac = eta
        ac = w :* g :* exp(etac :- JK.cen)
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
// when the model has tied failure times and was fitted with exact ties (tiemode 2), or with
// efron ties and weights other than 1 (tiemode 1; not allowed by stcox), so that Stata refits
// stcox instead (efron with tied failure times: Efron likelihood).
// nfname: number of clusters whose Newton fit failed; their codes go to the Stata matrix flname and the replicate PAFs of the others to the external
// __stpaf_PF (TAB unchanged), for Stata to refit those clusters and call _stpaf_jackfill
void _stpaf_jackfast(real matrix X, real colvector t0, real colvector t1, real colvector d,
    real colvector g, real colvector w, real colvector clus, real colvector strat,
    real rowvector b, real scalar nsc, real rowvector tau, real scalar level,
    real scalar logci, string scalar tabname, real scalar warm, real scalar tiemode,
    string scalar itname, string scalar nfname, string scalar flname, string scalar vpname)
{
    struct stpaf_dat scalar D
    struct stpaf_jk scalar JK
    pointer() scalar pPF
    pointer(real matrix) scalar pV, pM, pC
    real matrix PF, Q1, Q0, mt, DX, M
    real colvector cnt, wj, Lc, fl
    real rowvector bj, b0
    real scalar G, nt, s, t, gg, it, itsum, nok, rcg

    D = _stpaf_prep(X, t0, t1, d, g, w, clus, strat, b)
    JK = _stpaf_jkset(X, D)
    if (tiemode) {
        // Efron/exact likelihoods differ from Breslow only at tied failure times
        cnt = _stpaf_bsum(J(rows(JK.of), 1, 1), JK.inf, JK.kf, JK.Ktot)
        if (max(cnt) > 1) {
            // the Efron code assumes unit weights (0 for the left-out cluster)
            if (tiemode == 1 & all(w :== 1)) JK.efron = 1
            else {
                st_numscalar(itname, -1)
                st_numscalar(nfname, 0)
                return
            }
        }
    }
    // counterfactual columns: differences of the stored versions from the observed columns
    pV = findexternal("__stpaf_XV")
    pM = findexternal("__stpaf_XM")
    pC = findexternal("__stpaf_XC")
    if (pV == NULL | pM == NULL | pC == NULL) _error("stpaf internal error: counterfactual columns not found")
    M = *pM
    if (cols(*pC)) DX = (*pV) - X[., *pC]
    else DX = J(rows(X), 0, .)
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
        PF[gg, .] = _stpaf_jkpaf(JK, X, wj, g, bj, Lc, nsc, Q1, Q0, DX, M)
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
    _stpaf_jkvar(PF, level, logci, tabname, vpname)
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

// -------------------------------------------- average sequential (Shapley) PAFs (0.12.0)
// elements of s separated by blanks or tabs outside parentheses -> Stata locals pfx1, pfx2, ...;
// pfx + "n" = number of elements, pfx + "bad" = 1 if the parentheses are unbalanced
void _stpaf_toplocals(string scalar s, string scalar pfx)
{
    real scalar i, dep, n, bad
    string scalar ch, cur

    n = 0
    dep = 0
    bad = 0
    cur = ""
    for (i = 1; i <= strlen(s); i++) {
        ch = substr(s, i, 1)
        if (ch == "(") dep++
        else if (ch == ")") {
            dep--
            if (dep < 0) bad = 1
        }
        if ((ch == " " | ch == char(9)) & dep == 0) {
            if (cur != "") {
                n++
                st_local(pfx + strofreal(n), cur)
                cur = ""
            }
        }
        else cur = cur + ch
    }
    if (cur != "") {
        n++
        st_local(pfx + strofreal(n), cur)
    }
    if (dep != 0) bad = 1
    st_local(pfx + "n", strofreal(n))
    st_local(pfx + "bad", strofreal(bad))
}

// 1 if s is "(" ... ")" with the first parenthesis closed by the last character
real scalar _stpaf_isgroup(string scalar s)
{
    real scalar i, dep, n
    string scalar ch

    n = strlen(s)
    if (n < 2) return(0)
    if (substr(s, 1, 1) != "(" | substr(s, n, 1) != ")") return(0)
    dep = 0
    for (i = 1; i <= n; i++) {
        ch = substr(s, i, 1)
        if (ch == "(") dep++
        else if (ch == ")") {
            dep--
            if (dep == 0 & i < n) return(0)
        }
    }
    return(dep == 0)
}

// subsets of K factors (rows, 0/1 by factor): row 1 = all factors (joint scenario), then the
// other non-empty subsets by size and, within size, in lexicographic order of the factor
// numbers; Shapley weights W (K x subsets): factor k's average over the K! orders of the extra
// PAF from adding it, sum over S containing k of |S-1|!(K-|S|)!/K! [PAF(S) - PAF(S - k)], i.e.
// w_kS = (|S|-1)!(K-|S|)!/K! if k in S and -|S|!(K-|S|-1)!/K! otherwise (PAF(empty) = 0);
// sum_k w_kS = 1 for the full set and 0 for every other subset
void _stpaf_shapmat(real scalar K, string scalar subname, string scalar wname)
{
    real scalar ns, m, s, k, sz
    real matrix SUB, Z, W

    ns = 2^K - 1
    SUB = J(ns, K, 0)
    for (m = 1; m <= ns; m++) {
        for (k = 1; k <= K; k++) SUB[m, k] = mod(floor(m/2^(k-1)), 2)
    }
    Z = sort((rowsum(SUB), SUB), (1, -(2..(K+1))))
    SUB = Z[., 2..(K+1)]
    SUB = SUB[ns, .] \ SUB[1..(ns-1), .]
    W = J(K, ns, 0)
    for (s = 1; s <= ns; s++) {
        sz = sum(SUB[s, .])
        for (k = 1; k <= K; k++) {
            if (SUB[s, k]) W[k, s] = factorial(sz - 1)*factorial(K - sz)/factorial(K)
            else W[k, s] = -factorial(sz)*factorial(K - sz - 1)/factorial(K)
        }
    }
    st_matrix(subname, SUB)
    st_matrix(wname, W)
}

// final table of combined(): T (subset scenarios x horizons, joint first) and their covariance
// V (influence function or jackknife) -> rows: joint, the K Shapley PAFs, the other subsets.
// A Shapley row is W[k, .] applied to the subset PAFs at the same horizon (as are its
// influence function and jackknife replicates), so its variance is a V a'. CI on the PAF scale,
// normal quantile (dfq missing) or t(dfq); E = observed expected events, E_cf and se_logR
// missing. V is replaced by A V A' for the rows of the final table.
void _stpaf_shapfin(string scalar tabname, string scalar vpname, string scalar wname,
    real scalar nt, real scalar level, real scalar dfq)
{
    real matrix T, V, W, A, O
    real scalar K, ns, nin, nout, k, s, t, r, q, b, se

    T = st_matrix(tabname)
    V = st_matrix(vpname)
    W = st_matrix(wname)
    K = rows(W)
    ns = cols(W)
    nin = ns*nt
    if (rows(T) != nin | rows(V) != nin | cols(V) != nin) {
        _error("stpaf internal error: subset scenarios do not match the Shapley weights")
    }
    nout = (ns + K)*nt
    A = J(nout, nin, 0)
    O = J(nout, cols(T), .)
    q = (dfq < . ? invttail(dfq, (1 - level/100)/2) : invnormal(1 - (1 - level/100)/2))
    for (t = 1; t <= nt; t++) {
        A[t, t] = 1
        O[t, .] = T[t, .]
        for (s = 2; s <= ns; s++) {
            r = (K + s - 1)*nt + t
            A[r, (s-1)*nt + t] = 1
            O[r, .] = T[(s-1)*nt + t, .]
            O[r, 1] = K + s
        }
        for (k = 1; k <= K; k++) {
            r = k*nt + t
            for (s = 1; s <= ns; s++) A[r, (s-1)*nt + t] = W[k, s]
        }
    }
    V = A*V*A'
    V = (V + V')/2
    for (k = 1; k <= K; k++) {
        for (t = 1; t <= nt; t++) {
            r = k*nt + t
            b = A[r, .]*T[., 3]
            se = sqrt(max((V[r, r], 0)))
            O[r, .] = (k + 1, T[t, 2], b, se, b - q*se, b + q*se, T[t, 7], ., ., T[t, 10])
        }
    }
    st_matrix(tabname, O)
    st_matrix(vpname, V)
}

// ------------------------------------------------------------ multiple imputation
// completed-data command of an mi estimate command line: the text after the first colon
// that is outside parentheses and double quotes (saving(C:\...) has its colon inside)
string scalar _stpaf_micmd(string scalar cl)
{
    real scalar i, dep, inq
    string scalar ch

    dep = 0
    inq = 0
    for (i = 1; i <= strlen(cl); i++) {
        ch = substr(cl, i, 1)
        if (ch == char(34)) inq = !inq
        else if (!inq) {
            if (ch == "(") dep++
            else if (ch == ")") dep--
            else if (ch == ":" & dep == 0) return(strtrim(substr(cl, i + 1, .)))
        }
    }
    return("")
}

// file of the per-imputation fits named in an mi estimate command line: the file of the
// saving() option (abbreviated to sav at least), else the file after using, in the mi part of
// the line (before the first colon outside parentheses and quotes); "" if none.
// Checked in Python (mifile_parse_check.py, line-by-line port).
string scalar _stpaf_mifile(string scalar cl)
{
    real scalar n, i, j, k, dep, inq, op, cp
    string scalar ch, pre, msk, L, c, f
    string rowvector forms

    // mi part of the command line
    n = strlen(cl)
    pre = cl
    dep = 0
    inq = 0
    for (i = 1; i <= n; i++) {
        ch = substr(cl, i, 1)
        if (ch == char(34)) inq = !inq
        else if (!inq) {
            if (ch == "(") dep++
            else if (ch == ")") dep--
            else if (ch == ":" & dep == 0) {
                pre = substr(cl, 1, i - 1)
                break
            }
        }
    }
    // mask everything inside quotes or parentheses
    n = strlen(pre)
    msk = ""
    dep = 0
    inq = 0
    for (i = 1; i <= n; i++) {
        ch = substr(pre, i, 1)
        if (ch == ")" & !inq) dep--
        msk = msk + ((dep == 0 & !inq) ? ch : "#")
        if (ch == char(34)) inq = !inq
        else if (ch == "(" & !inq) dep++
    }
    L = " " + subinstr(strlower(msk), ",", " ", .)
    // saving(...): position of its opening parenthesis in pre
    forms = ("saving(", "savin(", "savi(", "sav(")
    op = 0
    for (k = 1; k <= 4; k++) {
        j = strpos(L, " " + forms[k])
        if (j) {
            op = j + strlen(forms[k]) - 1
            break
        }
    }
    if (op) {
        dep = 0
        inq = 0
        cp = 0
        for (i = op; i <= n; i++) {
            ch = substr(pre, i, 1)
            if (ch == char(34)) inq = !inq
            else if (!inq) {
                if (ch == "(") dep++
                else if (ch == ")") {
                    dep--
                    if (dep == 0) {
                        cp = i
                        break
                    }
                }
            }
        }
        if (cp == 0) return("")
        c = substr(pre, op + 1, cp - op - 1)
        // file name: up to the first comma outside quotes (suboption replace)
        f = c
        inq = 0
        for (i = 1; i <= strlen(c); i++) {
            ch = substr(c, i, 1)
            if (ch == char(34)) inq = !inq
            else if (ch == "," & !inq) {
                f = substr(c, 1, i - 1)
                break
            }
        }
        f = strtrim(f)
        if (strlen(f) >= 2) {
            if (substr(f, 1, 1) == char(34) & substr(f, strlen(f), 1) == char(34)) {
                f = substr(f, 2, strlen(f) - 2)
            }
        }
        return(f)
    }
    // mi estimate using filename
    j = strpos(L, " using ")
    if (j == 0) return("")
    i = j + 6
    while (i <= n) {
        if (substr(pre, i, 1) != " ") break
        i++
    }
    if (i > n) return("")
    if (substr(pre, i, 1) == char(34)) {
        k = strpos(substr(pre, i + 1, .), char(34))
        if (k == 0) return("")
        return(substr(pre, i + 1, k - 1))
    }
    f = ""
    while (i <= n) {
        ch = substr(pre, i, 1)
        if (ch == " " | ch == "," | ch == ":") break
        f = f + ch
        i++
    }
    return(f)
}

// add .ster when the file name has no extension
string scalar _stpaf_sterfn(string scalar f)
{
    return(pathsuffix(f) == "" ? f + ".ster" : f)
}

// command line normalized for comparison (blanks)
string scalar _stpaf_cmdnorm(string scalar s)
{
    string scalar t

    t = stritrim(strtrim(s))
    t = subinstr(t, " ,", ",", .)
    t = subinstr(t, ", ", ",", .)
    t = subinstr(t, "( ", "(", .)
    t = subinstr(t, " )", ")", .)
    return(t)
}

// Rubin's rules for each column of Q (M x k completed-data estimates) with variances U;
// nuc: complete-data degrees of freedom (missing: large sample, Rubin 1987; otherwise
// Barnard and Rubin 1999, as mi estimate). Returns k x 6: estimate, total variance T,
// df (missing = infinite), relative variance increase, fraction of missing information,
// within-imputation variance.
real matrix _stpaf_rubin(real matrix Q, real matrix U, real scalar nuc)
{
    real scalar M, k, c, qb, ub, b, t, r, nul, nuo, nu, gam, fmi
    real matrix out

    M = rows(Q)
    k = cols(Q)
    if (M < 2) _error("at least 2 imputations are required")
    out = J(k, 6, .)
    for (c = 1; c <= k; c++) {
        qb = mean(Q[., c])
        ub = mean(U[., c])
        b = sum((Q[., c] :- qb):^2)/(M - 1)
        t = ub + (1 + 1/M)*b
        if (t <= 0) {
            // no variability at all (e.g. a scenario that does not change the model)
            r = 0
            nu = .
            fmi = 0
        }
        else if (ub <= 0) {
            // all variability between imputations
            r = .
            nu = M - 1
            fmi = 1
        }
        else {
            if (b <= 0) {
                r = 0
                nul = .
            }
            else {
                r = (1 + 1/M)*b/ub
                nul = (M - 1)*(1 + 1/r)^2
            }
            if (nuc < .) {
                gam = (1 + 1/M)*b/t
                nuo = nuc*(nuc + 1)*(1 - gam)/(nuc + 3)
                nu = (nul < . ? 1/(1/nul + 1/nuo) : nuo)
                fmi = 1 - ((nu + 1)/(nu + 3))/((nuc + 1)/(nuc + 3))*ub/t
            }
            else {
                nu = nul
                fmi = (nul < . ? (r + 2/(nul + 3))/(r + 1) : 0)
            }
        }
        out[c, .] = (qb, t, nu, r, fmi, ub)
    }
    return(out)
}

// combination of the per-imputation r(table)s stacked in allname (M blocks of k rows,
// columns scenario horizon paf se lb ub E E_cf se_logR events): Rubin's rules on log(1-PAF)
// (logci) or on PAF; the Shapley rows of combined() (scenarios 2..K+1) always on PAF.
// Output k x 13 (the 10 columns plus df, riv, fmi); pfname: M x k PAFs, sename: M x k standard
// errors on the combination scale of each row; vsname: sum over the imputations of the
// within-imputation covariances of the PAFs (k x k), vpname: total covariance on the PAF scale,
// mean within + (1 + 1/M) between
void _stpaf_micomb(string scalar allname, real scalar M, real scalar nuc, real scalar level,
    real scalar logci, real scalar K, string scalar vsname, string scalar outname,
    string scalar pfname, string scalar sename, string scalar vpname)
{
    real matrix A, P, SP, SL, E, E0, EV, RL, RP, out, SEM
    real rowvector Rc
    real colvector sc, frc, lg
    real scalar k, c, al, q, paf, se, selr, lb, ub, ul

    A = st_matrix(allname)
    if (mod(rows(A), M)) _error("stpaf internal error: per-imputation tables differ in size")
    k = rows(A)/M
    P  = colshape(A[., 3], k)
    SP = colshape(A[., 4], k)
    E  = colshape(A[., 7], k)
    E0 = colshape(A[., 8], k)
    SL = colshape(A[., 9], k)
    EV = colshape(A[., 10], k)
    if (hasmissing(P)) _error("stpaf internal error: PAF missing in an imputation")
    sc = A[1..k, 1]
    frc = (K > 0 ? ((sc :>= 2) :& (sc :<= K + 1)) : J(k, 1, 0))
    lg = selectindex(frc :== 0)
    RL = J(k, 6, .)
    if (rows(lg)) {
        if (max(P[., lg]) < 1) RL[lg, .] = _stpaf_rubin(ln(1 :- P[., lg]), SL[., lg]:^2, nuc)
        else if (logci) _error(459, "PAF >= 1 in an imputation: log(1-PAF) is undefined; use citrans(none)")
    }
    RP = _stpaf_rubin(P, SP:^2, nuc)
    SEM = (logci ? SL : SP)
    al = (1 - level/100)/2
    out = J(k, 13, .)
    out[., 1..2] = A[1..k, 1..2]
    for (c = 1; c <= k; c++) {
        ul = (logci & !frc[c])
        Rc = (ul ? RL[c, .] : RP[c, .])
        q = (Rc[3] < . ? invttail(Rc[3], al) : invnormal(1 - al))
        selr = (frc[c] ? . : sqrt(RL[c, 2]))
        if (ul) {
            paf = 1 - exp(Rc[1])
            se = (1 - paf)*selr
            lb = 1 - exp(Rc[1] + q*selr)
            ub = 1 - exp(Rc[1] - q*selr)
        }
        else {
            paf = Rc[1]
            se = sqrt(Rc[2])
            lb = paf - q*se
            ub = paf + q*se
        }
        if (frc[c]) SEM[., c] = SP[., c]
        out[c, 3..13] = (paf, se, lb, ub, mean(E[., c]),
            (hasmissing(E0[., c]) ? . : mean(E0[., c])), selr,
            mean(EV[., c]), Rc[3], Rc[4], Rc[5])
    }
    st_matrix(outname, out)
    st_matrix(pfname, P)
    st_matrix(sename, SEM)
    st_matrix(vpname, st_matrix(vsname)/M + (1 + 1/M)*variance(P))
}

end
