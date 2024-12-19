{smcl}
{* *! version 1.0.0 19dec2024}{...}
{vieweralsosee "[ST] disco" "help disco"}{...}
{vieweralsosee "[ST] disco_estat" "help disco_estat"}{...}
{viewerjumpto "Syntax" "disco_plot##syntax"}{...}
{viewerjumpto "Description" "disco_plot##description"}{...}
{viewerjumpto "Options" "disco_plot##options"}{...}
{viewerjumpto "Examples" "disco_plot##examples"}{...}

{title:Title}

{phang}
{bf:disco_plot} {hline 2} Generate plots after DiSCo (Distributional Synthetic Controls)

{marker description}
{title:Description}

{pstd}
{cmd:disco_plot} creates visualizations of distributional treatment effects after {cmd:disco} 
estimation. It can display quantile functions, CDFs, and their differences over time, with 
optional confidence intervals.

{marker syntax}
{title:Syntax}

{p 8 17 2}
{cmdab:disco_plot}{cmd:,} {opt agg(string)} {opt m(integer)} {opt g(integer)} {opt t_max(integer)} {opt doci(integer)} {opt cl(real)}

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt agg(string)}}type of plot ("quantile", "quantileDiff", "cdf", or "cdfDiff"){p_end}
{synopt:{opt m(integer)}}number of quantile points{p_end}
{synopt:{opt g(integer)}}number of grid points for CDF{p_end}
{synopt:{opt t_max(integer)}}maximum time period{p_end}
{synopt:{opt doci(integer)}}whether to display confidence intervals (0/1){p_end}
{synopt:{opt cl(real)}}confidence level (e.g., 95){p_end}
{synoptline}
{p2colreset}{...}

{marker options}
{title:Options}

{dlgtab:Main}

{phang}
{opt agg(string)} specifies the type of plot to generate. Options are:

{p2colset 9 28 30 2}
{p2col:"quantile"}plot treated vs synthetic quantile functions{p_end}
{p2col:"quantileDiff"}plot differences in quantiles{p_end}
{p2col:"cdf"}plot treated vs synthetic CDFs{p_end}
{p2col:"cdfDiff"}plot differences in CDFs{p_end}

{phang}
{opt m(integer)} specifies the number of quantile points to use in the plot.

{phang}
{opt g(integer)} specifies the number of grid points to use for CDF plots.

{phang}
{opt t_max(integer)} specifies the maximum time period to plot.

{phang}
{opt doci(integer)} specifies whether to display confidence intervals (1) or not (0).

{phang}
{opt cl(real)} specifies the confidence level as a percentage (e.g., 95 for 95% confidence intervals).

{marker examples}
{title:Examples}

{pstd}Plot quantile differences with confidence intervals:{p_end}
{phang2}{cmd:. disco_plot, agg(quantileDiff) m(100) g(100) t_max(10) doci(1) cl(95)}{p_end}

{pstd}Plot CDFs without confidence intervals:{p_end}
{phang2}{cmd:. disco_plot, agg(cdf) m(100) g(100) t_max(10) doci(0) cl(95)}{p_end}

{pstd}Plot quantile functions:{p_end}
{phang2}{cmd:. disco_plot, agg(quantile) m(100) g(100) t_max(10) doci(0) cl(95)}{p_end}

{marker author}
{title:Author}

{pstd}
David Van Dijcke{break}
University of Michigan, Ann Arbor{break}
{browse "mailto:dvdijcke@umich.edu":dvdijcke@umich.edu}