# stpaf
A Stata command that estimates the population attributable fraction (PAF) of possibly time-varying covariates after `stcox`, for single events (Cox) and recurrent events (Lin-Wei-Yang-Ying), including Ghosh-Lin data prepared by `ghoshlin`. Confidence intervals come from an influence function or an exact cluster jackknife (the default with fewer than 100 clusters).

To install the command, run the following command in Stata:
##
    net from "https://szabobarnaki.github.io/stpaf/"

## Example
    stset tstop, id(id) failure(event) time0(tstart) exit(time .)
    stcox i.smoke age bmi, vce(cluster id)
    stpaf, eliminate(smoke) horizon(1 2 3)

See `help stpaf` for the options, methods and stored results. `stpaf_cert.do` (installed with `net get stpaf`) runs the certification tests.
