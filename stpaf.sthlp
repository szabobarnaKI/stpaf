{smcl}
{* *! version 0.4.0  02nov2026}{...}
{vieweralsodisplay}
{viewerjumpto "Syntax" "stpaf##syntax"}{...}
{viewerjumpto "Description" "stpaf##description"}{...}
{viewerjumpto "Options" "stpaf##options"}{...}
{viewerjumpto "Methods" "stpaf##methods"}{...}
{viewerjumpto "Examples" "stpaf##examples"}{...}
{viewerjumpto "Stored results" "stpaf##results"}{...}
{viewerjumpto "References" "stpaf##references"}{...}
{title:Title}

{phang}
{bf:stpaf} {hline 2} Population attributable fraction after {helpb stcox}
(Cox and Lin-Wei-Yang-Ying models), allowing time-varying covariates,
with influence-function or cluster-jackknife confidence intervals


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:stpaf} [{it:if}] [{it:in}] {cmd:,}
{{opt elim:inate(varlist)} | {opt set(var=# [var=# ...])} | {opt cap(var=# [var=# ...])}}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Counterfactual scenario (at least one required)}
{synopt:{opt elim:inate(varlist)}}set each variable to 0 (e.g. eliminate an exposure){p_end}
{synopt:{opt set(var=# ...)}}set each variable to the given value; no spaces around {cmd:=}{p_end}
{synopt:{opt cap(var=# ...)}}replace {it:var} by min({it:var}, #) (theoretical-minimum-risk level for continuous exposures){p_end}
{synopt:{opt each}}in addition to the joint scenario, report each specified change on its own{p_end}

{syntab:Reporting}
{synopt:{opt hor:izon(numlist)}}horizon(s) {it:tau} in analysis time; default is the whole follow-up{p_end}
{synopt:{opt l:evel(#)}}confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt citr:ans(log|none)}}scale used for the CI; default {cmd:log}, i.e. log(1-PAF){p_end}

{syntab:Variance}
{synopt:{opt vce(if|jackknife)}}influence function (default) or exact cluster jackknife{p_end}
{synopt:{opt cl:uster(varname)}}cluster variable; default {cmd:e(clustvar)}, then the {cmd:stset} {opt id()}{p_end}
{synopt:{opt nowarm}}{cmd:vce(jackknife)}: start each re-estimation from zero instead of {cmd:e(b)}{p_end}
{synopt:{opt refit}}{cmd:vce(jackknife)}: refit {cmd:stcox} for each cluster (slower; method of version 0.3){p_end}
{synopt:{opt nodots}}suppress the replication dots of {cmd:vce(jackknife)} with {opt refit}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
The data must be {helpb stset}. {cmd:stpaf} is a postestimation command for {helpb stcox}.
{it:if}/{it:in} restrict the person-time over which the attributable fraction is
evaluated, see {help stpaf##methods:Methods}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:stpaf} estimates the proportion of events that would have been avoided
had one or more (possibly time-varying) covariates been changed to a counterfactual
value over the whole follow-up, using the model fitted by {cmd:stcox}.
Both the Cox model (one event per subject) and the marginal-mean model of
Lin, Wei, Yang and Ying (2000) for recurrent events (a counting-process {cmd:stset}
with {opt id()} and {cmd:stcox, vce(cluster }{it:id}{cmd:)}) are covered, because the two
share the same estimating equations; only the interpretation changes
(hazard-weighted events versus expected number of recurrent events).

{pstd}
Confidence intervals come from a closed-form influence function that accounts for the
uncertainty in the regression coefficients, in the Breslow baseline hazard and
in the covariate distribution, and that sums contributions within clusters.
No resampling is needed; the cost is a few passes through the data
(O(n log n) for the sweeps).

{pstd}
Time-varying covariates are used exactly as they are stored in the (start, stop] records
of the {cmd:stset} data.  A counterfactual that changes a covariate changes it at all times.
Derived variables (cumulative exposure, lags, interaction with time created with
{helpb stsplit}) are not recomputed: set every variable that carries the exposure.


{marker options}{...}
{title:Options}

{phang}
{opt eliminate(varlist)}, {opt set()}, {opt cap()} define the counterfactual scenario.  They
are applied to the data in memory (temporarily) and the design matrix of the {cmd:stcox}
model is rebuilt, so factor-variable terms and interactions are handled correctly.
All changes given together form the {it:joint} scenario.

{phang}
{opt each} additionally reports the PAF for each change alone (partial PAFs;
these do not add up to the joint PAF in general).

{phang}
{opt horizon(numlist)} reports the cumulative PAF over (0, {it:tau}] for each {it:tau}.
Regression coefficients are always those of the full model.

{phang}
{opt citrans(log|none)}: with {cmd:log} the interval is
1 - (1-PAF)exp(±z se[log(1-PAF)]), which keeps the limits below 1.

{phang}
{opt vce(if|jackknife)} chooses the variance estimator.  {cmd:vce(if)}, the default, uses the
closed-form influence function (instantaneous).  {cmd:vce(jackknife)} re-estimates the model
G times, each time leaving out one cluster (same terms, strata and {cmd:stset} weights),
recomputes the PAF on each reduced sample, and uses the jackknife variance
(G-1)/G sum (theta_(-g) - mean)^2 of log(1-PAF) for the log-scale CI and of PAF for the
reported standard error.  It is recommended when there are fewer than about 100 clusters
(see {help stpaf##performance:Performance}).
{cmd:vce(jackknife)} is not available with {cmd:fweight}s and no clustering.

{pmore}
By default the reduced-sample estimates are computed inside {cmd:stpaf} by Newton-Raphson
for the Breslow partial likelihood, giving the left-out cluster weight zero, so that
{cmd:stcox} is not refitted.  Each fit starts from {cmd:e(b)} and is iterated until the
maximum relative change in the coefficients and the scaled gradient U'I^(-1)U are both below
1e-12.  This criterion is much tighter than the default convergence criteria of {cmd:stcox};
it is needed because the starting values are already close to the solution, so that a relative
change in the log likelihood (Stata's {cmd:ltolerance()}) can be met after a single step that
has not yet reached the optimum.  With it, the solution does not depend on the starting values
(see {opt nowarm}), and the estimates are the same as those of {cmd:stcox} refits up to the
convergence tolerance of {cmd:stcox}.  The method is used when the model was fitted with
{cmd:breslow} ties (the default of {cmd:stcox}) or when no failure times are tied (then
{cmd:efron} and exact ties give the same likelihood); otherwise {cmd:stcox} is refitted
and a note is displayed.  The number of clusters, the mean number of iterations and the
method are shown in the output and stored in {cmd:r()}.

{phang}
{opt nowarm} starts each jackknife re-estimation from zero instead of from {cmd:e(b)}.  The
results are the same (to about 1e-12); the option is provided to verify this and costs about
50% more iterations.

{phang}
{opt refit} refits {cmd:stcox} for each left-out cluster instead (the method of {cmd:stpaf}
0.3), using {cmd:stcox}'s own convergence criteria.  It is slower (G {cmd:stcox} fits) and is
provided for checking.

{phang}
{opt nodots} suppresses the replication dots displayed by {cmd:vce(jackknife)} {opt refit}.

{phang}
{opt cluster(varname)} defines the clusters (independent units) for both variance
estimators; {cmd:vce(jackknife)} leaves out one of these clusters at a time.  It overrides
{cmd:e(clustvar)} and the {cmd:stset} {opt id()}.  Because {cmd:stpaf} computes the variance
itself from the information matrix and the cluster scores, only the cluster definition of the
{cmd:stcox} fit matters, not its {cmd:vce()}: {cmd:stcox}, {cmd:vce(robust)} followed by
{cmd:stpaf}, {cmd:cluster(}{it:id}{cmd:)} gives the same results as {cmd:stcox}, {cmd:vce(cluster} {it:id}{cmd:)}
followed by {cmd:stpaf}.  Point estimates do not depend on the clusters.  Use it when the clusters
are not declared to {cmd:stset}, in particular after {cmd:ghoshlin} (see
{help stpaf##weights:Weights}), or for clusters above the subject level (e.g. centres).
Subjects should be nested within the clusters; a note is displayed when the {cmd:stset}
{opt id()} or {cmd:e(clustvar)} is not nested within {cmd:cluster()}.  Missing values of the
cluster variable in the estimation sample are not allowed.


{marker methods}{...}
{title:Methods and formulas}

{pstd}
Let {it:N_i(t)} count the events of cluster {it:i}, {it:Y_i(t)} be the at-risk indicator,
{it:X_i(t)} the covariates and {it:X_i^c(t)} their counterfactual values.  With
the Breslow estimator {it:dL(t) = dN.(t)/S0(t)}, {it:S0(t) = sum_j Y_j exp(b'X_j)},
the expected numbers of events are

{p 8 8 2}
E = sum_i integral Y_i exp(b'X_i) dL, {space 4} E0 = sum_i integral Y_i exp(b'X_i^c) dL,

{pstd}
and PAF(tau) = 1 - E0(tau)/E(tau).  With no {it:if}, E equals the observed number of events.
With an {it:if} condition, the sums over {it:i} are restricted to the person-time
records satisfying the condition (for example the exposed person-time, giving an
attributable fraction among the exposed), whereas the baseline hazard still uses all records.

{pstd}
Using the martingale residuals {it:dM_i = dN_i - Y_i exp(b'X_i) dL}, the influence function
of log(E0/E) for cluster {it:i} is

{p 8 8 2}
phi_i = psi0_i/E0 - psi_i/E, {space 4} psi_i = integral h dM_i + integral T_i dL - E/G + D'(I^-1 U_i),

{pstd}
where {it:h = T/S0}, {it:T_i} the cluster's contribution to the (restricted) risk-set sum,
{it:U_i} its partial-likelihood score, {it:I} the information matrix and
{it:D = dE/db} at fixed event times (psi0 is analogous with the counterfactual
covariates).  The variance estimate is G/(G-1) times the sum of phi_i^2 over the G clusters.
Clusters are defined by {opt cluster()}, otherwise by {cmd:e(clustvar)} from {cmd:stcox},
otherwise by the {cmd:stset} {opt id()}, otherwise by record.

{pstd}
{bf:Caveats.}  (1) The estimand is the hazard/mean-function based attributable fraction.
For Cox models with competing or terminal events it is not a risk-based PAF.  (2) Valid
causal reading requires the usual no-unmeasured-confounding and correct model assumptions;
the CI covers sampling variability only.  (3) The variance uses Breslow-type score and
information at {cmd:e(b)}; with {cmd:efron} or {cmd:exactp} ties and many tied
event times it is an approximation (a note is displayed when the score at {cmd:e(b)} is not zero).
(4) {cmd:tvc()}, {cmd:offset()} and {cmd:shared()} are not supported.


{marker performance}{...}
{title:Performance of the confidence intervals}

{pstd}
In simulations of recurrent events with a gamma frailty, two strata, a time-varying
exposure and a correctly specified LWYY model (1000 replications for the default, 500 for the
jackknife; Monte Carlo SE of coverage 0.007-0.010), the point estimate was unbiased at all
sample sizes, and the coverage of nominal 95% log-scale intervals was:

{p2colset 8 30 32 2}{...}
{p2col:Clusters}influence function / jackknife{p_end}
{p2col:50}0.90 / 0.93{p_end}
{p2col:100}0.94 / 0.95{p_end}
{p2col:400}0.94 / -{p_end}
{p2col:1000}0.95-0.96 / -{p_end}
{p2colreset}{...}

{pstd}
With 50 clusters the clustered (sandwich) SE of the regression coefficient reported by
{cmd:stcox} itself covered 0.89, i.e. the influence-function CI of {cmd:stpaf} behaves like
the model's own robust variance.  The jackknife SE matched the empirical SD already with
50 clusters.  {cmd:stpaf} displays a note when the default is used with fewer than 100 clusters.

{marker weights}{...}
{title:Weights}

{pstd}
Weights are taken from {helpb stset} ({cmd:fweight}, {cmd:pweight}, {cmd:iweight}) and may differ
between the records of a subject.  They enter the Breslow estimator, the expected
numbers of events and the influence function as known, record-level multipliers.
Consequently, a weighted fit that uses time-varying weights
(for example inverse-probability-of-censoring weights created by a data-preparation command and
supplied with {cmd:stset} {it:[pw=weight]}) is handled; the weights are treated as fixed, so uncertainty from
estimating them is {it:not} reflected in the standard errors.

{pstd}
Data prepared by {cmd:ghoshlin} are {cmd:stset} {it:[pw=weight_c]} without {opt id()}, because {cmd:stset}
does not allow probability weights that vary within subjects.  The subject identifier must then
be given either to {cmd:stcox} ({cmd:vce(cluster} {it:id}{cmd:)}) or to {cmd:stpaf} ({cmd:cluster(}{it:id}{cmd:)});
after frequency weights expanded by {cmd:ghoshlin}, use {cmd:repid} instead of {it:id}.  If neither is given,
{cmd:stpaf} stops with an error rather than treat the records of a subject as independent.

{pstd}
With {cmd:fweight}s and no clustering (no {opt cluster()}, no {cmd:vce(cluster)} and no {cmd:stset} {opt id()}) a record with weight
{it:w} is treated as {it:w} independent replicates; the result equals that from the expanded data set.
In all other cases (including {cmd:fweight}s with clusters) the weighted cluster sums are used, as for {cmd:pweight}s.
Weights must not be negative.

{marker examples}{...}
{title:Examples}

{pstd}Cox model, baseline covariates{p_end}
{phang2}{cmd:. webuse drugtr, clear}{p_end}
{phang2}{cmd:. stcox i.drug age, breslow}{p_end}
{phang2}{cmd:. stpaf, set(drug=0)}{p_end}

{pstd}Recurrent events (LWYY) with a time-varying exposure; {cmd:exit(time .)} keeps follow-up after the first event{p_end}
{phang2}{cmd:. stset tstop, id(id) failure(event) time0(tstart) exit(time .)}{p_end}
{phang2}{cmd:. stcox i.smoke age bmi, vce(cluster id)}{p_end}
{phang2}{cmd:. stpaf, eliminate(smoke) horizon(1 2 3) each}{p_end}

{pstd}Few clusters: exact cluster jackknife{p_end}
{phang2}{cmd:. stpaf, eliminate(smoke) vce(jackknife)}{p_end}

{pstd}Ghosh-Lin model (recurrent events with a terminal event) after {cmd:ghoshlin}; run {cmd:stpaf} in the transformed frame{p_end}
{phang2}{cmd:. sysuse ghoshlin_example, clear}{p_end}
{phang2}{cmd:. stset fup_days, id(ID) failure(events == 1, 2, 3) exit(time .)}{p_end}
{phang2}{cmd:. ghoshlin, failure("inlist(events,1,3)") competing("inlist(events,1,2)") keep(age sex diabetes) frame(transformed,replace)}{p_end}
{phang2}{cmd:. frame change transformed}{p_end}
{phang2}{cmd:. stcox c.age i.sex i.diabetes, vce(cluster ID) breslow}{p_end}
{phang2}{cmd:. stpaf, eliminate(diabetes) horizon(365 730)}{p_end}

{pstd}Same, with the subject identifier given to {cmd:stpaf} instead of {cmd:stcox}{p_end}
{phang2}{cmd:. stcox c.age i.sex i.diabetes, breslow}{p_end}
{phang2}{cmd:. stpaf, eliminate(diabetes) horizon(365 730) cluster(ID)}{p_end}

{pstd}Continuous exposure reduced to a threshold, in the exposed person-time only{p_end}
{phang2}{cmd:. stpaf if sbp > 140, cap(sbp=140)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}{cmd:stpaf} stores in {cmd:r()} results of the joint scenario at the last horizon:

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(paf)}}population attributable fraction{p_end}
{synopt:{cmd:r(se)}}standard error of PAF{p_end}
{synopt:{cmd:r(lb)}, {cmd:r(ub)}}confidence limits{p_end}
{synopt:{cmd:r(E)}, {cmd:r(E_cf)}}expected events, observed and counterfactual{p_end}
{synopt:{cmd:r(se_logR)}}standard error of log(1-PAF){p_end}
{synopt:{cmd:r(N_clust)}}number of clusters{p_end}
{synopt:{cmd:r(N_fail)}}number of failures{p_end}
{synopt:{cmd:r(jk_iter)}}{cmd:vce(jackknife)} (Newton): mean number of iterations per cluster{p_end}
{p2col 5 16 20 2: Macros}{p_end}
{synopt:{cmd:r(cfactual)}}counterfactual scenario{p_end}
{synopt:{cmd:r(citrans)}}CI transformation{p_end}
{synopt:{cmd:r(vce)}}{cmd:if} or {cmd:jackknife}{p_end}
{synopt:{cmd:r(clustvar)}}cluster variable (empty if clustered on records){p_end}
{synopt:{cmd:r(jk_method)}}{cmd:vce(jackknife)}: {cmd:newton}, {cmd:newton, nowarm} or {cmd:refit}{p_end}
{p2col 5 16 20 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}one row per scenario and horizon; columns scenario, horizon, paf, se, lb, ub, E, E_cf, se_logR, events{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}Lin, D. Y., L. J. Wei, I. Yang, and Z. Ying. 2000. Semiparametric regression
for the mean and rate functions of recurrent events. {it:Journal of the Royal Statistical Society, Series B} 62: 711-730.{p_end}

{phang}Chen, L., D. Y. Lin, and D. Zeng. 2010. Attributable fraction functions for
censored event times. {it:Biometrika} 97: 713-726.{p_end}
