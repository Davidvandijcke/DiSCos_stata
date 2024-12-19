{smcl}
{* *! version 1.0.0 19dec2024}{...}
{vieweralsosee "[ST] synth_runner" "help synth_runner"}{...}
{viewerjumpto "Syntax" "disco##syntax"}{...}
{viewerjumpto "Description" "disco##description"}{...}
{viewerjumpto "Options" "disco##options"}{...}
{viewerjumpto "Examples" "disco##examples"}{...}
{viewerjumpto "Method" "disco##method"}{...}
{viewerjumpto "Stored results" "disco##results"}{...}
{viewerjumpto "References" "disco##references"}{...}
{viewerjumpto "Author" "disco##author"}{...}

{title:Title}

{phang}
{bf:disco} {hline 2} Distributional Synthetic Controls.

{marker description}
{title:Description}

{pstd}
{cmd:disco} implements the Distributional Synthetic Controls (DiSCo) method based on 
Gunsilius (2023), extending the synthetic 
control methodology of Abadie and Gardeazabal (2003) and Abadie, Diamond, and Hainmueller (2010) 
to distributions. Instead of focusing solely on aggregate mean outcomes, DiSCo constructs 
synthetic control weights that replicate an entire outcome distribution of a treated unit 
from a set of control units.

{pstd}
By reproducing the quantile function (or, optionally, mixture of distributions) of the 
treated unit before treatment, DiSCo identifies a set of weights that can be used to 
form a synthetic control distribution in all time periods. This synthetic distribution 
then serves as an estimate of the counterfactual distribution of the treated unit 
in the absence of treatment. By comparing the observed treated distribution with 
this synthetic counterpart, DiSCo estimates distributional treatment effects such as 
differences in quantiles, CDFs, and other distributional functionals.

{pstd}
This approach is suitable for settings where one aims to identify heterogeneous treatment 
effects along the entire distribution of outcomes, rather than focusing solely on averages. 
The method is particularly useful when richer micro-level data are available at the unit level 
(e.g., states, firms) but individuals within those units cannot be tracked over time.

{pstd}
{cmd:disco} also supports bootstrap inference for confidence intervals, permutation tests 
analogous to the classical synthetic control permutation inference, and graphical summaries 
of results.

{pstd}
Please cite Gunsilius (2023) and Van Dijcke, Gunsilius, and Wright (2024) when using this package.

{marker syntax}
{title:Syntax}

{p 8 17 2}
{cmd:disco} {it:varlist(3)} [{it:if}] [{it:in}], {opt idtarget(#)} {opt t0(#)} [{it:options}]

{pstd}
{it:varlist} must contain exactly three variables in the following order:

{phang2}1. Outcome variable (numeric){p_end}
{phang2}2. Unit ID variable (numeric){p_end}
{phang2}3. Time period variable (integer){p_end}

{pstd}
{cmd:idtarget()} and {cmd:t0()} are required.

{marker options}
{title:Options}

{dlgtab:Required}

{phang}
{opt idtarget(#)} specify the id of the treated unit.

{phang}
{opt t0(#)} specify the first treatment period.

{dlgtab:Optional}

{phang}
{opt m(integer)} number of quantile points used for quantile-based estimation. default is 100.

{phang}
{opt g(integer)} number of grid points for cdf estimation. default is 100.

{phang}
{opt ci} compute bootstrap confidence intervals for distributional effects.

{phang}
{opt boots(integer)} number of bootstrap replications for confidence intervals. default is 300.

{phang}
{opt cl(real)} confidence level for intervals. default is 0.95.

{phang}
{opt qmin(real)} minimum quantile for estimation range. default is 0.

{phang}
{opt qmax(real)} maximum quantile for estimation range. default is 1.

{phang}
{opt nosimplex} do not constrain weights to lie in a unit simplex. by default, weights are nonnegative 
and sum to one. specifying {cmd:nosimplex} allows weights to take any values that sum to one.

{phang}
{opt mixture} use the mixture (cdf-based) approach instead of the quantile-based approach.

{phang}
{opt permutation} perform a permutation test by treating each control unit as a "placebo" treated unit 
and computing test statistics. returns a p-value.

{phang}
{opt seed(integer)} set the random seed for reproducibility. default is -1 (no seed set).

{phang}
{opt nouniform} when computing confidence intervals, do not compute uniform confidence bands; 
only pointwise intervals are computed.

{phang}
{opt agg(string)} specify the type of aggregation for summary statistics. one of:
{p_end}
{phang2}- {cmd:"quantile"}: summarize estimated quantile functions{p_end}
{phang2}- {cmd:"cdf"}: summarize estimated cdfs{p_end}
{phang2}- {cmd:"quantilediff"}: summarize differences in quantiles between treated and synthetic{p_end}
{phang2}- {cmd:"cdfdiff"}: summarize differences in cdfs between treated and synthetic{p_end}

{phang}
{opt samples(numlist)} specify quantile or cdf points for summary statistics. for quantiles, these are in [0,1]. 
for cdfs, these are values of the outcome variable.

{phang}
{opt graph} produce graphical output of the results by time period.

{marker examples}
{title:Examples}

{pstd}Basic usage with confidence intervals:{p_end}
{phang2}{cmd:. disco y id time, idtarget(1) t0(3) m(50) g(100) ci boots(200) cl(0.90)}{p_end}

{pstd}Using mixture approach:{p_end}
{phang2}{cmd:. disco outcome unit t, idtarget(2) t0(10) mixture ci}{p_end}

{pstd}With permutation test and graphing:{p_end}
{phang2}{cmd:. disco wage county year, idtarget(10) t0(2005) permutation seed(12345) graph}{p_end}

{marker method}
{title:Method and Formulation}

{pstd}
Distributional synthetic controls extend the idea of synthetic controls to the entire 
outcome distribution. Instead of matching average outcomes, we match entire quantile 
functions or CDFs of control units to replicate the pre-treatment distribution of a 
treated unit. Post-treatment differences then yield distributional treatment effects.

{dlgtab:Quantile-based (2-Wasserstein) approach}

{pstd}
Consider a treated unit indexed by 1 and control units indexed by j=2,...,J+1 observed 
over periods t=1,...,T, with T0 < T as the last pre-treatment period. Let Y_{j t} be 
outcomes for unit j in period t. We want to estimate the counterfactual distribution 
Y_{1 t, N} that would have prevailed for the treated unit in the absence of treatment.

{pstd}
The key object is the quantile function of Y_{j t}, denoted F_{Y_{j t}}^{-1}(q). By 
forming a weighted average of these quantile functions for the control units:

{phang}F_{Y_{1 t, N}}^{-1}(q) = ∑_{j=2}^{J+1} λ_j^* F_{Y_{j t}}^{-1}(q){p_end}

{pstd}
we obtain a synthetic control distribution for the treated unit. The optimal weights 
λ_j^* are chosen to minimize the 2-Wasserstein distance between the treated 
unit's pre-treatment quantile functions and a weighted combination of the controls' 
pre-treatment quantile functions:

{phang}λ_{t}^* = arg min_{λ ∈ Δ^{J}} ∫_{0}^{1} |F_{Y_{1 t}}^{-1}(q) - ∑_{j=2}^{J+1}λ_j F_{Y_{j t}}^{-1}(q)|^2 dq.{p_end}

{dlgtab:CDF-based (1-Wasserstein) approach}

{pstd}
In some cases, it may be more natural to replicate the treated unit's distribution function 
F_{Y_{1 t}}(y) itself rather than its quantile function. Instead of solving for weights 
that match quantiles, one can find weights to match cumulative distribution functions (CDFs):

{phang}λ_{t}^* = arg min_{λ ∈ Δ^{J}} ∫_{ℝ} |F_{Y_{1 t}}(y) - ∑_{j=2}^{J+1} λ_j F_{Y_{j t}}(y)| dy.{p_end}

{pstd}
This corresponds to using the 1-Wasserstein distance (or equivalently L^1-distance between 
distributions). The 1-Wasserstein approach mixes entire distributions directly.

{marker results}
{title:Stored results}

{pstd}
{cmd:disco} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(t0)}}first treatment period{p_end}
{synopt:{cmd:e(cl)}}confidence level used{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}"disco"{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(agg)}}aggregation type, if specified{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(weights)}}estimated synthetic control weights{p_end}
{synopt:{cmd:e(quantile_diff)}}differences in quantiles by time{p_end}
{synopt:{cmd:e(cdf_diff)}}differences in cdfs by time{p_end}
{synopt:{cmd:e(quantile_synth)}}synthetic quantiles{p_end}
{synopt:{cmd:e(quantile_t)}}treated unit quantiles{p_end}
{synopt:{cmd:e(cdf_synth)}}synthetic cdfs{p_end}
{synopt:{cmd:e(cdf_t)}}treated unit cdfs{p_end}
{synopt:{cmd:e(summary_stats)}}summary statistics if agg option specified{p_end}

{pstd}If {cmd:ci} specified:{p_end}
{synopt:{cmd:e(qdiff_lower)}}lower ci bounds for quantile differences{p_end}
{synopt:{cmd:e(qdiff_upper)}}upper ci bounds for quantile differences{p_end}
{synopt:{cmd:e(cdiff_lower)}}lower ci bounds for cdf differences{p_end}
{synopt:{cmd:e(cdiff_upper)}}upper ci bounds for cdf differences{p_end}

{marker related}
{title:Additional Commands}

{phang}{cmd:disco_estat}: summarize aggregated statistics if specified with agg() option.{p_end}

{phang}{cmd:disco_plot}: generate plots for quantiles or cdfs across time.{p_end}

{marker references}
{title:References}

{phang}
Abadie, Alberto, and Javier Gardeazabal. 2003. "The Economic Costs of Conflict: A Case Study of the Basque Country."
{browse "http://dx.doi.org/10.1257/000282803321455188":American Economic Review 93(1): 113–132.}
{p_end}

{phang}
Abadie, A. 2021. "Using Synthetic Controls: Feasibility, Data Requirements, and Methodological Aspects."
{browse "http://dx.doi.org/10.1257/jel.20191450":Journal of Economic Literature 59(2): 391-425.}
{p_end}

{phang}
Abadie, A., Diamond, A., & Hainmueller, J. 2010. "Synthetic Control Methods for Comparative Case Studies: Estimating the Effect of California's Tobacco Control Program."
{browse "http://dx.doi.org/10.1198/jasa.2009.ap08746":Journal of the American Statistical Association 105(490): 493-505.}
{p_end}

{phang}
Gunsilius, F. 2023. "Distributional Synthetic Controls."
{browse "http://dx.doi.org/10.3982/ECTA18260":Econometrica 91(3): 1105-1117.}
{p_end}

{phang}
Van Dijcke, D., Gunsilius, F., & Wright, A. L. 2024. "Return to Office and the Tenure Distribution."
{browse "https://bfi.uchicago.edu/working-paper/2024-56/":University of Chicago, Becker Friedman Institute for Economics Working Paper, (2024-56).}
{p_end}

{marker author}
{title:Author}

{pstd}
David Van Dijcke{break}
University of Michigan, Ann Arbor{break}
{browse "mailto:dvdijcke@umich.edu":dvdijcke@umich.edu}
{p_end}

{title:Version}

{pstd}
1.0.0 (December 2024)
{p_end}
