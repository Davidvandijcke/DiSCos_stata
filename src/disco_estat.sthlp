{smcl}
{* *! version 1.0.0 19dec2024}{...}
{vieweralsosee "[ST] disco" "help disco"}{...}
{vieweralsosee "[ST] disco_plot" "help disco_plot"}{...}
{viewerjumpto "Syntax" "disco_estat##syntax"}{...}
{viewerjumpto "Description" "disco_estat##description"}{...}
{viewerjumpto "Examples" "disco_estat##examples"}{...}
{viewerjumpto "Stored results" "disco_estat##results"}{...}

{title:Title}

{phang}
{bf:disco_estat} {hline 2} Post-estimation statistics for DiSCo (Distributional Synthetic Controls)

{marker description}
{title:Description}

{pstd}
{cmd:disco_estat} displays summary statistics after {cmd:disco} estimation. It provides a detailed 
summary of distributional treatment effects, including point estimates, standard errors, and 
confidence intervals across different time periods.

{marker syntax}
{title:Syntax}

{p 8 17 2}
{cmdab:disco_estat} {cmd:summary}

{pstd}
{cmd:disco_estat} requires that {cmd:disco} has been run previously with the {cmd:agg()} option.

{marker examples}
{title:Examples}

{pstd}Display summary statistics after running disco with quantile differences:{p_end}
{phang2}{cmd:. disco y id time, idtarget(1) t0(3) agg(quantileDiff)}{p_end}
{phang2}{cmd:. disco_estat summary}{p_end}

{pstd}Display summary statistics for CDF differences:{p_end}
{phang2}{cmd:. disco y id time, idtarget(1) t0(3) agg(cdfDiff)}{p_end}
{phang2}{cmd:. disco_estat summary}{p_end}

{marker results}
{title:Stored results}

{pstd}
{cmd:disco_estat} requires the following to be stored in {cmd:e()} from a previous {cmd:disco} estimation:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}"disco"{p_end}
{synopt:{cmd:e(agg)}}type of aggregation ("quantileDiff" or "cdfDiff"){p_end}
{synopt:{cmd:e(cl)}}confidence level{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(summary_stats)}}matrix containing time period, range, effect, std. err., and confidence intervals{p_end}

{marker author}
{title:Author}

{pstd}
David Van Dijcke{break}
University of Michigan, Ann Arbor{break}
{browse "mailto:dvdijcke@umich.edu":dvdijcke@umich.edu}