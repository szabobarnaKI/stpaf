{smcl}
{* *! version 0.14.0  04nov2026}{...}
{viewerjumpto "Syntax" "stpaf##syntax"}{...}
{viewerjumpto "Description" "stpaf##description"}{...}
{viewerjumpto "Options" "stpaf##options"}{...}
{viewerjumpto "Methods" "stpaf##methods"}{...}
{viewerjumpto "Level-specific PAFs" "stpaf##levels"}{...}
{viewerjumpto "Average sequential (Shapley) PAFs" "stpaf##shapley"}{...}
{viewerjumpto "Multiple imputation" "stpaf##mi"}{...}
{viewerjumpto "Performance" "stpaf##performance"}{...}
{viewerjumpto "Weights" "stpaf##weights"}{...}
{viewerjumpto "Examples" "stpaf##examples"}{...}
{viewerjumpto "Stored results" "stpaf##results"}{...}
{viewerjumpto "References" "stpaf##references"}{...}
{viewerjumpto "AI disclosure" "stpaf##ai"}{...}
{title:Title}

{p2colset 5 15 17 2}{...}
{p2col :{cmd:stpaf} {hline 2}}Population attributable fraction after {helpb stcox} (Cox and Lin-Wei-Yang-Ying models), allowing time-varying covariates, with influence-function or cluster-jackknife confidence intervals, also after {helpb mi estimate}{cmd:: stcox}{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:stpaf} [{it:if}] [{it:in}] {cmd:,}
{{opt elim:inate(varlist)} | {opt set(var=# [var=# ...])} | {opt cap(var=# [var=# ...])}}
[{opt each} {it:options}]

{p 8 16 2}
{cmd:stpaf} [{it:if}] [{it:in}] {cmd:,}
{opt varlev:els(varname[=#])}
[{it:options}]

{p 8 16 2}
{cmd:stpaf} [{it:if}] [{it:in}] {cmd:,}
{opt comb:ined(factorlist)} [{opt sub:sets} {it:options}]

{p 4 4 2}
where {it:factorlist} is two to eight factors, each

{p 8 16 2}
{it:element} | {cmd:(}{it:element} [{it:element} ...]{cmd:)}

{p 4 4 2}
and {it:element} is {opt elim:inate(varlist)}, {opt set(var=# [var=# ...])} or
{opt cap(var=# [var=# ...])}.  Each variable of an element outside parentheses is one factor;
the elements inside one pair of parentheses together form one factor.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Counterfactual scenario (one of the three forms above is required)}
{synopt:{opt elim:inate(varlist)}}set each variable to 0 (e.g. eliminate an exposure){p_end}
{synopt:{opt set(var=# ...)}}set each variable to the given value; no spaces around {cmd:=}{p_end}
{synopt:{opt cap(var=# ...)}}replace {it:var} by min({it:var}, #) (theoretical-minimum-risk level for continuous exposures){p_end}
{synopt:{opt each}}in addition to the joint scenario, report each specified change on its own{p_end}
{synopt:{opt varlev:els(varname[=#])}}joint PAF of moving all levels of the factor variable {it:varname} to the reference level ({it:#}; default the base level), and the PAF of each level moved alone{p_end}
{synopt:{opt comb:ined(factorlist)}}joint PAF of 2 to 8 factors and the average sequential (Shapley) PAF of each factor{p_end}
{synopt:{opt sub:sets}}with {opt combined()}: also display the PAFs of the other subsets of factors{p_end}

{syntab:Reporting}
{synopt:{opt hor:izon(numlist)}}horizon(s) {it:tau} in analysis time; default is the whole follow-up{p_end}
{synopt:{opt l:evel(#)}}confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt citr:ans(log|none)}}scale of the confidence interval; default {cmd:log}, i.e. log(1-PAF){p_end}

{syntab:Variance}
{synopt:{opt vce(if|jackknife)}}influence function or exact cluster jackknife; default {cmd:jackknife} with fewer than 100 clusters, otherwise {cmd:if}{p_end}
{synopt:{opt cl:uster(varname)}}cluster variable; default {cmd:e(clustvar)}, then the {cmd:stset} {opt id()}{p_end}
{synopt:{opt nowarm}}jackknife: start each re-estimation from zero instead of {cmd:e(b)}{p_end}
{synopt:{opt refit}}jackknife: refit {cmd:stcox} for each cluster (slower; for checking){p_end}
{synopt:{opt nodots}}suppress the replication dots of {opt refit} and the imputation dots after {cmd:mi estimate}{p_end}

{syntab:After mi estimate}
{synopt:{opt mifile(filename)}}file of the {cmd:stcox} fits saved by {cmd:mi estimate, saving()}; default: the file named in {cmd:e(cmdline_mi)}{p_end}
{synopt:{opt mirefit}}refit {cmd:stcox} in each imputation even when saved fits are available{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
The data must be {helpb stset} (or {helpb mi stset} for imputed data). {cmd:stpaf} is a postestimation
command for {helpb stcox} and for {helpb mi estimate}{cmd:: stcox}.
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
(O(n log n) for the sweeps).  With fewer than 100 clusters, where the influence-function
interval is too narrow, the default is an exact cluster jackknife instead
(see {opt vce()}).

{pstd}
Besides the PAF of one (joint) change, {cmd:stpaf} reports partial PAFs of each change
({opt each}), level-specific PAFs of a categorical exposure that add up to its joint PAF
({opt varlevels()}), and average sequential (Shapley) PAFs of several risk factors that add up
to their joint PAF ({opt combined()}).

{pstd}
Time-varying covariates are used exactly as they are stored in the (start, stop] records
of the {cmd:stset} data.  A counterfactual that changes a covariate changes it at all times.
Derived variables (cumulative exposure, lags, interaction with time created with
{helpb stsplit}) are not recomputed: set every variable that carries the exposure.

{pstd}
After {cmd:mi estimate: stcox}, {cmd:stpaf} combines the estimates of the imputed data sets
by Rubin's rules; see {help stpaf##mi:Multiple imputation}.


{marker options}{...}
{title:Options}

{dlgtab:Counterfactual scenario}

{phang}
{opt eliminate(varlist)}, {opt set()}, {opt cap()} define the counterfactual scenario.  They
are applied to the data in memory (temporarily) and the design matrix of the {cmd:stcox}
model is rebuilt, so factor-variable terms and interactions (also with continuous variables,
e.g. {cmd:c.z#c.age}) are handled correctly.  All changes given together form the
{it:joint} scenario.

{phang}
{opt each} additionally reports the PAF for each change alone (partial PAFs;
these do not add up to the joint PAF in general; see {opt combined()} for a decomposition that does).

{phang}
{opt varlevels(varname[=#])} reports level-specific PAFs of a categorical exposure that enters
the model as a factor variable ({cmd:i.}{it:varname}, in main effects and/or interactions).  The
first scenario (labelled {it:varname}{cmd:=}{it:#}) moves every record to the reference level {it:#}
(the same as {cmd:set(}{it:varname}{cmd:=}{it:#}{cmd:)}); then, for each other level {it:k}, the
scenario labelled {it:k}{cmd:.}{it:varname} moves only the records at level {it:k} to the reference,
leaving all other records unchanged.  As for the other scenarios, the design matrix is rebuilt, so
interactions are recomputed.  The level PAFs add up exactly to the joint PAF; see
{help stpaf##levels:Level-specific PAFs}.  The reference level defaults to the base level of
{it:varname} in {cmd:e(b)} (e.g. 3 after {cmd:ib3.}{it:varname}); {it:varname}{cmd:=}{it:#} sets it to any
level of {it:varname} in {cmd:e(b)}, and is required when the base level cannot be determined (e.g.
{cmd:ibn.}).  The levels are those of {it:varname} in the column names of {cmd:e(b)}, not those found in
the data, so that they are the same in every imputation after {cmd:mi estimate}; a level without
records in the {it:if} subset has PAF 0.  {opt varlevels()} cannot be combined with {opt eliminate()},
{opt set()}, {opt cap()}, {opt each} or {opt combined()}.

{phang}
{opt combined(factorlist)} reports the joint PAF of changing K risk factors (2 <= K <= 8) and, for
each factor, its average sequential PAF (Eide and Gefeller 1995), i.e. its Shapley value: the extra
PAF from also changing that factor, averaged over all K! orders in which the factors can be changed.
The K Shapley PAFs add up exactly to the joint PAF.  A factor is defined by the elements
{opt eliminate()}, {opt set()} and {opt cap()} (with the same meaning as the options of that name):
each variable of an element outside parentheses is one factor, and the elements inside one pair of
parentheses together form one factor, e.g.

{pmore2}
{cmd:combined((eliminate(d1) set(d2=1)) cap(bmi=25) set(sbp=120))}

{pmore}
has three factors, labelled {cmd:d1+d2}, {cmd:bmi} and {cmd:sbp}, and {cmd:combined(eliminate(smoke diab))}
has two, {cmd:smoke} and {cmd:diab}.  Groups cannot be nested; a variable may appear only once; blanks
around {cmd:=} are allowed; {cmd:eliminate} can be abbreviated to {cmd:elim}.  {cmd:stpaf} computes the
PAFs of all 2^K - 1 non-empty subsets of the factors (each a scenario with its design matrix rebuilt,
so interactions within and between factors are recomputed), and each Shapley PAF is a fixed linear
combination of them.  See {help stpaf##shapley:Average sequential (Shapley) PAFs}.  {opt combined()}
cannot be combined with {opt eliminate()}, {opt set()}, {opt cap()}, {opt each} or {opt varlevels()}.

{phang}
{opt subsets} (with {opt combined()} only) also displays the PAFs of the 2^K - 2 subsets of factors
other than the joint scenario.  They are always stored in {cmd:r(table)} and {cmd:r(V_paf)}.

{dlgtab:Reporting}

{phang}
{opt horizon(numlist)} reports the cumulative PAF over (0, {it:tau}] for each {it:tau}.
Regression coefficients are always those of the full model.

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence intervals.
The default is {cmd:level(95)} or as set by {helpb set level}.

{phang}
{opt citrans(log|none)} sets the scale of the confidence interval.  With {cmd:log} (the default)
the interval is 1 - (1-PAF)exp(±q se[log(1-PAF)]), which keeps the limits below 1; with {cmd:none}
it is PAF ± q se(PAF) (q is the normal or t quantile, see {opt vce()}).  After {cmd:mi estimate},
{opt citrans()} also sets the scale on which Rubin's rules are applied.  The Shapley PAFs of
{opt combined()} always use the PAF scale.

{dlgtab:Variance}

{phang}
{opt vce(if|jackknife)} chooses the variance estimator.  {cmd:vce(if)} uses the
closed-form influence function (no re-estimation).  {cmd:vce(jackknife)} re-estimates the model
G times, each time leaving out one cluster (same terms, strata and {cmd:stset} weights),
recomputes the PAF on each reduced sample, and uses the jackknife variance
(G-1)/G sum (theta_(-g) - mean)^2 of log(1-PAF) for the log-scale CI and of PAF for the
reported standard error.  The jackknife confidence limits use the t distribution with G-1
degrees of freedom (as Stata's {helpb jackknife} prefix does); the influence-function limits use
the normal distribution.
{cmd:vce(jackknife)} is not available with {cmd:fweight}s and no clustering.

{pmore}
If {opt vce()} is not specified, the jackknife is used when there are fewer than 100 clusters
(where the influence-function CI is too narrow; see {help stpaf##performance:Performance})
and the influence function otherwise.  Either method can be requested explicitly whatever the
number of clusters; with {cmd:vce(if)} and fewer than 100 clusters a note is displayed.  When
the jackknife was chosen by default and cannot be computed (for example, an {it:if} restriction
that lies within a single cluster leaves a reduced sample without expected events), {cmd:stpaf}
displays a note with the error code and reports the influence-function results; an explicitly
requested {cmd:vce(jackknife)} stops with the error instead.

{pmore}
By default the reduced-sample estimates are computed inside {cmd:stpaf} by Newton-Raphson
for the partial likelihood of the model (Breslow, or Efron for {cmd:efron} ties), giving the
left-out cluster weight zero, so that
{cmd:stcox} is not refitted.  Each fit starts from {cmd:e(b)} and is iterated until the
maximum relative change in the coefficients and the scaled gradient U'I^(-1)U are both below
1e-12.  This criterion is much tighter than the default convergence criteria of {cmd:stcox};
it is needed because the starting values are already close to the solution, so that a relative
change in the log likelihood (Stata's {cmd:ltolerance()}) can be met after a single step that
has not yet reached the optimum.  With it, the solution does not depend on the starting values
(see {opt nowarm}), and the estimates are the same as those of {cmd:stcox} refits up to the
convergence tolerance of {cmd:stcox}.  The method is used when the model was fitted with
{cmd:breslow} ties (the default of {cmd:stcox}) or {cmd:efron} ties, and with {cmd:exactm} or
{cmd:exactp} ties when no failure times are tied (the likelihoods are then the same as
Breslow's); with exact ties and tied failure times {cmd:stcox} is refitted for every cluster and
a note is displayed.  As in the point estimate, the baseline hazard of each reduced sample is the
Breslow estimator at its coefficients, also after {cmd:efron}.  ({cmd:stcox} allows no {cmd:stset}
weights with {cmd:efron} ties, so weighted models always use the Newton jackknife.)  If the Newton
fit fails for a cluster (singular information matrix, e.g. a covariate that is constant once the
cluster is left out, or no convergence in 100 iterations), {cmd:stcox} is refitted for that cluster
only, and a note is displayed; {cmd:stcox} then omits collinear terms as usual.  The number of
clusters, the mean number of iterations, the number of {cmd:stcox} refits and the method are shown
in the output and stored in {cmd:r()}.

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

{phang}
{opt nowarm} starts each jackknife re-estimation from zero instead of from {cmd:e(b)}.  The
results are the same (to about 1e-12); the option is provided to verify this and costs about
50% more iterations.

{phang}
{opt refit} refits {cmd:stcox} for each left-out cluster instead, using {cmd:stcox}'s own
convergence criteria.  It is slower (G {cmd:stcox} fits) and is provided for checking.

{pmore}
{opt nowarm} and {opt refit} cannot be combined.  Without {opt vce()} either of them requests
the jackknife whatever the number of clusters; with {cmd:vce(if)} they are not allowed.

{phang}
{opt nodots} suppresses the replication dots displayed by {cmd:vce(jackknife)} {opt refit}
(and by the {cmd:stcox} refits of clusters whose Newton fit failed) and, after {cmd:mi estimate},
the dots displayed as the imputations are processed.

{dlgtab:After mi estimate}

{phang}
{opt mifile(filename)} names the file ({cmd:.ster}) with the completed-data {cmd:stcox} fits saved
by {cmd:mi estimate, saving(}{it:filename}{cmd:)}.  By default the file named in the {cmd:saving()}
option or after {cmd:using} in {cmd:e(cmdline_mi)} is used if it exists; specify {opt mifile()} when
the working directory has changed since {cmd:mi estimate}, or when the fits were saved by an earlier
{cmd:mi estimate} of the same model.  A file given in {opt mifile()} must exist.  See
{help stpaf##mi:Multiple imputation}.

{phang}
{opt mirefit} refits {cmd:stcox} in each imputation even when saved fits are available.  The
results are the same.


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
The influence function of PAF itself is -(1-PAF)phi_i.  The covariance of the PAFs of all rows of
{cmd:r(table)} (scenarios and horizons), stored in {cmd:r(V_paf)}, is G/(G-1) times the sum over the
clusters of the products of these influence functions (with the cluster jackknife: (G-1)/G times the
sum of the products of the centred replicate PAFs).  Its diagonal is the square of the reported
standard errors.

{pstd}
The term I^-1 U_i is the influence function of the coefficients.  For {cmd:breslow} ties, and for
any ties without tied failure times, U_i and I are those of the Breslow partial likelihood.  For
{cmd:efron} ties with tied failure times, U_i is the cluster sum of the
Efron score residuals (Therneau's centred residuals, as {cmd:predict, esr} after {cmd:stcox})
and I the Efron information, so that G/(G-1) sum (I^-1 U_i)(I^-1 U_i)' equals the
{cmd:vce(cluster)} variance of {cmd:stcox} (stored in {cmd:r(V_b)}); the other terms are exact
because the baseline hazard is the Breslow estimator at the Efron coefficients.  (In simulations
with many tied failure times, Breslow-type expressions would understate the variance of the
coefficients by about 10%.)

{pstd}
{bf:Computation.}  A column of the design matrix changes in a counterfactual scenario only if its
term contains a changed variable.  {cmd:stpaf} therefore stores only such columns, one version for
each distinct set of changes among the variables of the term, and builds all versions in as few
passes over the data as possible (one when the changed variables enter as main effects only).  The
risk-set sums at the event times are computed for all records at once, and the influence-function
covariance is accumulated over blocks of clusters, so that memory does not grow with the number of
clusters times the number of rows of {cmd:r(table)}.

{pstd}
{bf:Caveats.}  (1) The estimand is the hazard/mean-function based attributable fraction.
For Cox models with competing or terminal events it is not a risk-based PAF.  (2) Valid
causal reading requires the usual no-unmeasured-confounding and correct model assumptions;
the CI covers sampling variability only.  (3) With {cmd:exactm} or {cmd:exactp} ties and tied
failure times the influence function uses Breslow-type score and information at {cmd:e(b)}, an
approximation (a note is displayed when the score at {cmd:e(b)} is not zero); {cmd:vce(jackknife)}
refits the exact likelihood.  (4) {cmd:tvc()}, {cmd:offset()} and {cmd:shared()} are not supported.


{marker levels}{...}
{title:Level-specific PAFs}

{pstd}
Let x have levels r (reference), k_1, ..., k_L.  Moving the records at level k to r changes the
counterfactual risk score exp(b'X^c) of these records only; the baseline hazard dL is that of the
observed data in every scenario.  Hence E - E0_k is the sum, over the records at level k, of the
integral of {c -(}exp(b'X) - exp(b'X^c){c )-}Y dL, and summing over k gives E - E0 of the joint
scenario.  So PAF_k = (E - E0_k)/E, and PAF_k1 + ... + PAF_kL = PAF(joint) exactly, at every horizon,
with {it:if} restrictions, interactions and time-varying x (each record has one level; a subject
whose level changes contributes to several level PAFs).  The decomposition also holds in every
jackknife replicate, and for the influence function on the PAF scale.  {cmd:r(sumdiff)} is the
largest |sum of the level PAFs - joint PAF| over the horizons.

{pstd}
Standard errors and confidence intervals of each level PAF are computed as for any other scenario
(influence function or cluster jackknife, on the log(1-PAF) scale with {cmd:citrans(log)}).
Moving level k to the reference is the same as {cmd:eliminate(}{it:dk}{cmd:)} in an equivalent
model with indicator variables {it:dk} = (x == k) (and their interactions), which can be used as a check.

{pstd}
After {cmd:mi estimate}, Rubin's rules are applied to each row.  With {cmd:citrans(none)} (PAF scale)
the combined level PAFs add up exactly to the combined joint PAF; with {cmd:citrans(log)} (the default)
the combination is not linear in PAF and they add up only approximately (the difference is of the
order of the between-imputation variance); the difference is displayed and stored in {cmd:r(sumdiff)}.


{marker shapley}{...}
{title:Average sequential (Shapley) PAFs}

{pstd}
Partial PAFs of several risk factors ({opt each}) do not add up to the joint PAF, and the sequential
PAF of a factor (the extra PAF from changing it after some others) depends on the order of the changes.
The average sequential PAF of Eide and Gefeller (1995) averages the sequential PAF of factor k over
all K! orders; it is the Shapley value of factor k in the cooperative game with value PAF(S) for each
subset S of the factors:

{p 8 8 2}
PAF_k = sum over S not containing k of |S|!(K-|S|-1)!/K! {c -(}PAF(S + k) - PAF(S){c )-}, {space 2} PAF(empty) = 0.

{pstd}
Hence PAF_1 + ... + PAF_K = PAF(all factors) exactly (at every horizon and with {it:if} restrictions),
and PAF_k = sum_S w_kS PAF(S) with fixed weights w_kS (stored in {cmd:r(W_shap)}).  Each PAF(S) is
the PAF of the scenario that changes the factors in S, with the design matrix rebuilt, so
interactions between factors are handled (interactions are where the order matters).  Factors
are changed jointly at all times, as for any other scenario.

{pstd}
Because a Shapley PAF is a fixed linear combination of the subset PAFs, its influence function is
the same combination of their influence functions (on the PAF scale), and its jackknife replicates
the same combination of their replicates: Var(PAF_k) = w_k' V w_k with V the covariance of the
subset PAFs ({cmd:r(V_paf)}).  The decomposition holds exactly in every jackknife replicate.  As
PAF_k is a difference of PAFs and can be negative, its confidence interval is computed on the PAF
scale (normal quantile with the influence function, t(G-1) with the jackknife), whatever
{opt citrans()}; the joint and subset rows use {opt citrans()} as usual.  In {cmd:r(table)} the
Shapley rows have E (observed) but missing E_cf and se_logR; PAF_k x E is factor k's averaged share
of the prevented events.  {cmd:r(sumdiff)} is the largest |sum of the Shapley PAFs - joint PAF| over
the horizons.

{pstd}
The counterfactual design matrices of the 2^K - 1 subset scenarios share their columns (see
{help stpaf##methods:Computation}): a main effect of one factor needs one stored version, an
interaction of two factors three, so memory grows with the number of stored versions
({cmd:r(xc_cols)}, about N x 8 bytes each for N records), not with 2^K.  The computing time grows
roughly in proportion to 2^K.  For example, with K = 8 (255 scenarios) and 18 coefficients,
{cmd:stpaf} took about 7 s with the influence function and 13 s with the cluster jackknife
(20 clusters) on 20,000 records, and about 1 minute with the influence function on 125,000 records
(Stata/MP 19, 4 cores).

{pstd}
After {cmd:mi estimate}, the Shapley PAFs are computed in each imputation and combined by Rubin's
rules on the PAF scale (point estimate = mean over the imputations).  With {cmd:citrans(none)} they
add up exactly to the combined joint PAF; with {cmd:citrans(log)} the joint PAF is combined on the
log(1-PAF) scale, so they add up only approximately ({cmd:r(sumdiff)}, displayed).


{marker mi}{...}
{title:Multiple imputation}

{pstd}
After {cmd:mi estimate: stcox} (with or without {cmd:mi estimate}'s {opt post} option), {cmd:stpaf}
runs, for each imputation used by {cmd:mi estimate} ({cmd:e(m_est_mi)}), on the imputed data
({helpb mi extract}) with the options given (scenarios, horizons, {it:if}, clusters, variance
method).  Counterfactual changes of imputed covariates are therefore applied in each imputation.
{cmd:mi estimate}'s transformations, {opt nimputations()} and {opt imputations()} are respected
through {cmd:e(m_est_mi)}.

{pstd}
{bf:Completed-data fits.}  When {cmd:mi estimate} was run with {cmd:saving(}{it:miestfile}{cmd:)}
(or as {cmd:mi estimate using} {it:miestfile}, or the file is given in {opt mifile()}), the
{cmd:stcox} fits saved there are used, so that {cmd:stcox} is not refitted.  A saved fit is used for
an imputation only if (1) it is a fit of the same {cmd:stcox} command ({cmd:e(cmdline)}, blanks
ignored), (2) the estimation sample rebuilt from that command (the command's {it:if}/{it:in},
{cmd:_st}, no missing covariates, {opt strata()} or {cmd:vce(cluster)} variable) has {cmd:e(N)}
observations, and (3) its coefficients maximize the partial likelihood of the imputed data in
memory (largest standardized score at {cmd:e(b)} below 1e-4, computed by {cmd:stpaf}).  Otherwise,
and when no file is available, {cmd:stcox} is refitted in that imputation (the command is taken
from {cmd:e(cmdline_mi)}), and a note gives the reason.  Both ways give the same results.  Fits with
{cmd:exactm} or {cmd:exactp} ties and tied failure times cannot be checked by (3) and are always
refitted.  For imputation {it:m}, the {it:k}th of {cmd:e(m_est_mi)}, the fits numbered {it:k} and {it:m}
in the file are tried (the {it:k}th first when the file holds exactly as many fits as imputations used);
check (3) rejects a fit of another imputation.  The mean of the coefficients used must equal
{cmd:e(b_mi)}; otherwise the data or the imputations have changed since {cmd:mi estimate}
and {cmd:stpaf} stops (error 498).  Only the {cmd:stcox} fits are saved by {cmd:mi estimate}, so
the saving is one {cmd:stcox} fit per imputation; the within-imputation variance (influence
function or jackknife) is computed in every imputation.

{pstd}
The M estimates are combined by Rubin's rules on the scale of {opt citrans()}: with
{cmd:citrans(log)} (the default), Q_m = log(1-PAF_m) with within-imputation variance
se_logR_m^2; the reported PAF is 1 - exp(mean Q_m), the total variance is
T = W + (1+1/M)B (W the mean within-imputation variance, B the between-imputation variance),
{cmd:se_logR} = sqrt(T), the reported standard error of PAF is (1-PAF)sqrt(T), and the limits are
1 - exp(mean Q_m ± t sqrt(T)).  With {cmd:citrans(none)} the same rules are applied to PAF itself.
Shapley PAFs ({opt combined()}) are always combined on the PAF scale.

{pstd}
{cmd:r(V_paf)} after {cmd:mi estimate} is the multivariate version on the PAF scale: the mean of the
within-imputation covariances plus (1+1/M) times the between-imputation covariance of the PAFs.  Its
diagonal equals the squared standard errors of rows combined on the PAF scale; with {cmd:citrans(log)}
it differs slightly from them for the other rows (which are combined on the log(1-PAF) scale).

{pstd}
The t quantile uses degrees of freedom that depend on the within-imputation variance method,
as {cmd:mi estimate} does depending on whether the command stores residual degrees of freedom.
With the influence function (a large-sample method), the large-sample df of Rubin (1987),
(M-1)(1+1/r)^2 with r = (1+1/M)B/W, are used (normal quantile when B = 0).  With the cluster
jackknife, whose complete-data intervals use t(G-1), the small-sample df of Barnard and Rubin
(1999) are used with complete-data df G-1 (the smallest G if it varies between imputations).
Note that with B = 0 the Barnard-Rubin df are (G-1)G/(G+2), slightly below G-1, as in
{cmd:mi estimate}.  The relative variance increase r and the fraction of missing information
(the formulas of {cmd:mi estimate}) are reported as well.

{pstd}
The default variance method (jackknife with fewer than 100 clusters) is decided in each
imputation; if the default jackknife fails in some imputations but not in others,
{cmd:stpaf} stops and asks for an explicit {cmd:vce()}.  The jackknife is run in every
imputation, so its time is multiplied by M.  The weights, clusters and strata are those of each
imputed data set; the data must be {helpb mi stset}.

{pstd}
The validity of results obtained with {cmd:mi} data depends on the validity of the technique
used to handle the missing data.


{marker performance}{...}
{title:Performance of the confidence intervals}

{pstd}
In simulations of recurrent events with a gamma frailty, two strata, a time-varying
exposure and a correctly specified LWYY model (500-2000 replications; Monte Carlo SE of
coverage 0.005-0.010), the point estimate was essentially unbiased at all sample sizes
(relative bias at most 4% with 25 clusters).  The coverage of nominal 95% log-scale
intervals was (ranges: horizon 2 and end of follow-up, and repeated simulations; -: not studied):

{col 9}{hline 10}{c TT}{hline 44}
{col 9}{ralign 9:Clusters} {c |}{col 22}Influence{col 36}Jackknife,{col 50}Jackknife,
{col 9}{space 9} {c |}{col 22}function{col 36}normal{col 50}t(G-1)
{col 9}{hline 10}{c +}{hline 44}
{col 9}{ralign 9:25} {c |}{col 22}0.89{col 36}0.93-0.94{col 50}0.94-0.95
{col 9}{ralign 9:50} {c |}{col 22}0.90-0.92{col 36}0.93-0.95{col 50}0.94-0.95
{col 9}{ralign 9:100} {c |}{col 22}0.94{col 36}0.95{col 50}-
{col 9}{ralign 9:400} {c |}{col 22}0.94{col 36}-{col 50}-
{col 9}{ralign 9:1000} {c |}{col 22}0.95-0.96{col 36}-{col 50}-
{col 9}{hline 10}{c BT}{hline 44}

{pstd}
With 25 clusters the jackknife SE was about 6% larger than the empirical SD, but the misses
of its interval were one-sided (the interval lay below the true PAF), reflecting the skewness
of the estimator in small samples; the t(G-1) quantile recovers most of the shortfall.

{pstd}
With 50 clusters the clustered (sandwich) SE of the regression coefficient reported by
{cmd:stcox} itself covered 0.89, i.e. the influence-function CI of {cmd:stpaf} behaves like
the model's own robust variance, whereas the jackknife SE matched the empirical SD.  Therefore
the jackknife is the default with fewer than 100 clusters.

{pstd}
With multiply imputed covariates (Rubin's rules as described above, 5 or 20 imputations), the
coverage was 0.95-0.96 with 25 clusters and 0.94-0.95 with 50 clusters (default jackknife), and
0.94-0.98 with 400 clusters (influence function).

{pstd}
The Newton jackknife is fast: about 3 s for 400 clusters (1,641 recurrent events), and 12 s for 961
clusters and 35,002 records (the Ghosh-Lin example below), against 26 s and 452 s for
{cmd:stcox} refits.


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

{pstd}Exact cluster jackknife (the default with fewer than 100 clusters, e.g. clustering on centres){p_end}
{phang2}{cmd:. stpaf, eliminate(smoke) vce(jackknife)}{p_end}
{phang2}{cmd:. stpaf, eliminate(smoke) cluster(centre)}{p_end}

{pstd}Level-specific PAFs of a categorical, time-varying exposure (reference: the base level){p_end}
{phang2}{cmd:. stcox i.bmicat##c.age i.sex, vce(cluster id)}{p_end}
{phang2}{cmd:. stpaf, varlevels(bmicat) horizon(1 2 3)}{p_end}
{phang2}{cmd:. stpaf, varlevels(bmicat=2)}{p_end}

{pstd}Average sequential (Shapley) PAFs of three risk factors that interact; the two diabetes indicators form one factor{p_end}
{phang2}{cmd:. stcox i.smoke##c.sbp diab1 diab2 age, vce(cluster id)}{p_end}
{phang2}{cmd:. stpaf, combined(set(smoke=0) cap(sbp=120) (eliminate(diab1 diab2))) horizon(1 2 3)}{p_end}
{phang2}{cmd:. stpaf, combined(set(smoke=0) cap(sbp=120) (eliminate(diab1 diab2))) subsets}{p_end}
{phang2}{cmd:. matrix list r(V_paf)}{p_end}

{pstd}Continuous exposure reduced to a threshold, in the exposed person-time only{p_end}
{phang2}{cmd:. stpaf if sbp > 140, cap(sbp=140)}{p_end}

{pstd}Multiply imputed covariates (here bmi), recurrent events; the fits saved by {cmd:mi estimate} are reused{p_end}
{phang2}{cmd:. mi stset tstop, id(id) failure(event) time0(tstart) exit(time .)}{p_end}
{phang2}{cmd:. mi estimate, saving(miest, replace): stcox i.smoke age bmi, vce(cluster id)}{p_end}
{phang2}{cmd:. stpaf, eliminate(smoke) horizon(1 2 3)}{p_end}

{pstd}Later, from the saved fits{p_end}
{phang2}{cmd:. mi estimate using miest}{p_end}
{phang2}{cmd:. stpaf, eliminate(smoke) horizon(1 2 3)}{p_end}

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


{marker results}{...}
{title:Stored results}

{pstd}{cmd:stpaf} stores the following in {cmd:r()}; the scalars paf to se_logR refer to the joint
scenario (scenario 1) at the last horizon:

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(paf)}}population attributable fraction{p_end}
{synopt:{cmd:r(se)}}standard error of PAF{p_end}
{synopt:{cmd:r(lb)}, {cmd:r(ub)}}confidence limits{p_end}
{synopt:{cmd:r(E)}, {cmd:r(E_cf)}}expected events, observed and counterfactual{p_end}
{synopt:{cmd:r(se_logR)}}standard error of log(1-PAF){p_end}
{synopt:{cmd:r(N_clust)}}number of clusters{p_end}
{synopt:{cmd:r(N_fail)}}number of failures{p_end}
{synopt:{cmd:r(level)}}confidence level{p_end}
{synopt:{cmd:r(df)}}{cmd:vce(jackknife)}: degrees of freedom G-1 of the t quantile of the CI{p_end}
{synopt:{cmd:r(jk_iter)}}jackknife (Newton): mean number of iterations per cluster fitted by Newton{p_end}
{synopt:{cmd:r(jk_nrefit)}}jackknife (Newton): number of clusters refitted by {cmd:stcox} after a failed Newton fit{p_end}
{synopt:{cmd:r(reflevel)}}{opt varlevels()}: reference level{p_end}
{synopt:{cmd:r(K)}}{opt combined()}: number of factors{p_end}
{synopt:{cmd:r(sumdiff)}}{opt varlevels()}: largest |sum of the level PAFs - joint PAF| over the horizons; {opt combined()}: the same for the Shapley PAFs{p_end}
{synopt:{cmd:r(xc_cols)}}number of stored counterfactual column versions{p_end}
{synopt:{cmd:r(xc_builds)}}number of passes over the data to build them{p_end}
{synopt:{cmd:r(if_chunks)}}number of blocks of clusters over which the influence-function covariance was accumulated{p_end}

{p2col 5 16 20 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}{cmd:stpaf}{p_end}
{synopt:{cmd:r(cfactual)}}counterfactual scenario{p_end}
{synopt:{cmd:r(scen}{it:#}{cmd:)}}label of scenario {it:#} (rows of {cmd:r(table)} with scenario = {it:#}); with {opt combined()}: {cmd:joint}, {cmd:Shapley: }{it:factor}, or the factors of a subset separated by commas{p_end}
{synopt:{cmd:r(citrans)}}CI transformation{p_end}
{synopt:{cmd:r(vce)}}{cmd:if} or {cmd:jackknife}, the method actually used{p_end}
{synopt:{cmd:r(clustvar)}}cluster variable (empty if clustered on records){p_end}
{synopt:{cmd:r(jk_method)}}{cmd:vce(jackknife)}: {cmd:newton}, {cmd:newton, nowarm} or {cmd:refit}{p_end}
{synopt:{cmd:r(if_ties)}}score and information in the coefficient part of the influence function: {cmd:efron} (efron ties with tied failure times) or {cmd:breslow}{p_end}
{synopt:{cmd:r(levelvar)}}{opt varlevels()}: the factor variable{p_end}
{synopt:{cmd:r(levels)}}{opt varlevels()}: the levels moved to the reference, in the order of scenarios 2, 3, ...{p_end}
{synopt:{cmd:r(factors)}}{opt combined()}: factor labels (variables of a factor joined by {cmd:+}), in factor order{p_end}

{p2col 5 16 20 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}one row per scenario and horizon; columns scenario, horizon, paf, se, lb, ub, E, E_cf, se_logR, events.  With {opt combined()}: scenario 1 = joint, 2 to K+1 = Shapley PAFs of factors 1 to K, K+2 to K+2^K-1 = the other subsets (by size, then in factor order){p_end}
{synopt:{cmd:r(V_paf)}}covariance of the PAFs of all rows of {cmd:r(table)} (row names s{it:#}h{it:#}: scenario, horizon), from the influence function or the jackknife; its diagonal is se^2{p_end}
{synopt:{cmd:r(V_b)}}cluster-robust variance of the (non-omitted) coefficients used by the influence function, G/(G-1) sum (I^-1 U_i)(I^-1 U_i)'; equals {cmd:e(V)} of {cmd:stcox}, {cmd:vce(cluster)} on the same clusters (not stored with {cmd:fweight}s and no clustering){p_end}
{synopt:{cmd:r(W_shap)}}{opt combined()}: K x (2^K-1) weights of the Shapley PAFs; columns = the subset scenarios in {cmd:r(table)} order (scenario 1, then K+2, ...){p_end}
{synopt:{cmd:r(subsets)}}{opt combined()}: (2^K-1) x K 0/1 matrix, the factors in each subset scenario (rows as the columns of {cmd:r(W_shap)}){p_end}
{p2colreset}{...}

{pstd}After {cmd:mi estimate}, {cmd:r(E)}, {cmd:r(E_cf)} and {cmd:r(N_fail)} are means over the
imputations, {cmd:r(df)} holds the MI degrees of freedom (missing = large) and {cmd:r(jk_nrefit)} is
summed over the imputations; {cmd:r(V_b)}, {cmd:r(if_ties)}, {cmd:r(jk_iter)}, {cmd:r(xc_cols)},
{cmd:r(xc_builds)} and {cmd:r(if_chunks)} are not stored.  Additionally:

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(M)}}number of imputations used{p_end}
{synopt:{cmd:r(riv)}}relative variance increase due to missing data{p_end}
{synopt:{cmd:r(fmi)}}fraction of missing information{p_end}
{synopt:{cmd:r(df_c)}}jackknife: complete-data degrees of freedom (G-1){p_end}
{synopt:{cmd:r(mi_nsaved)}}number of imputations in which the saved {cmd:stcox} fit was used{p_end}

{p2col 5 16 20 2: Macros}{p_end}
{synopt:{cmd:r(mi)}}{cmd:mi}{p_end}
{synopt:{cmd:r(m_mi)}}imputation numbers used{p_end}
{synopt:{cmd:r(mi_fits)}}{cmd:saved} (all fits from the file), {cmd:refit} (none) or {cmd:mixed}{p_end}
{synopt:{cmd:r(mi_file)}}file of the saved fits (when it was found){p_end}

{p2col 5 16 20 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}as above, plus columns df, riv, fmi{p_end}
{synopt:{cmd:r(V_paf)}}mean within-imputation covariance + (1+1/M) between-imputation covariance of the PAFs (PAF scale){p_end}
{synopt:{cmd:r(paf_mi)}}PAF in each imputation (rows) for each row of {cmd:r(table)} (columns){p_end}
{synopt:{cmd:r(se_mi)}}within-imputation standard errors on the combination scale of each row (se_logR or se){p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}Barnard, J., and D. B. Rubin. 1999. Small-sample degrees of freedom with multiple imputation.
{it:Biometrika} 86: 948-955.{p_end}

{phang}Chen, L., D. Y. Lin, and D. Zeng. 2010. Attributable fraction functions for
censored event times. {it:Biometrika} 97: 713-726.{p_end}

{phang}Eide, G. E., and O. Gefeller. 1995. Sequential and average attributable fractions as aids in
the selection of preventive strategies. {it:Journal of Clinical Epidemiology} 48: 645-655.{p_end}

{phang}Lin, D. Y., L. J. Wei, I. Yang, and Z. Ying. 2000. Semiparametric regression
for the mean and rate functions of recurrent events. {it:Journal of the Royal Statistical Society, Series B} 62: 711-730.{p_end}

{phang}Rubin, D. B. 1987. {it:Multiple Imputation for Nonresponse in Surveys}. New York: Wiley.{p_end}

{phang}Therneau, T. M., and P. M. Grambsch. 2000. {it:Modeling Survival Data: Extending the Cox Model}. New York: Springer.{p_end}


{marker ai}{...}
{title:AI disclosure}

{pstd}
This command was developed using substantial work from Claude AI (Anthropic), with human supervision.
{p_end}

{marker author}{...}
{title:Author}

{pstd}
Barna Szabó-Söderberg{break}
{it:Karolinska Institutet}{break}
{ul:barna.szabo@ki.se}
{p_end}
