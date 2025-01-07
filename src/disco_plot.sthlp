{smcl}
{* *! version 1.0.0 19dec2024}{...}
{vieweralsosee "[ST] disco" "help disco"}{...}
{viewerjumpto "Syntax" "disco_plot##syntax"}{...}
{viewerjumpto "Description" "disco_plot##description"}{...}
{viewerjumpto "Options" "disco_plot##options"}{...}
{viewerjumpto "Examples" "disco_plot##examples"}{...}
{viewerjumpto "References" "disco_plot##references"}{...}
{viewerjumpto "Author" "disco_plot##author"}{...}

{title:disco_plot}

{phang}
{bf:disco_plot} {hline 2} Post-estimation plots for Distributional Synthetic Controls

{marker description}
{title:Description}

{pstd}
{cmd:disco_plot} creates visualizations after {cmd:disco} estimation. It can display quantile 
functions, CDFs, and their differences over time, with optional confidence intervals. The command 
automatically reads all necessary information from the stored results of the previous {cmd:disco} 
estimation.

{marker syntax}
{title:Syntax}

{p 8 17 2}
{cmd:disco_plot} [{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}

{marker options}
{title:Options}

{dlgtab:Plot Type}

{phang}
{opt agg(string)} type of plot to generate: (default: quantileDiff if not specified in disco)
{p_end}
{phang2}- {cmd:"quantile"}: plot treated vs synthetic quantile functions{p_end}
{phang2}- {cmd:"quantileDiff"}: plot differences in quantiles{p_end}
{phang2}- {cmd:"cdf"}: plot treated vs synthetic CDFs{p_end}
{phang2}- {cmd:"cdfDiff"}: plot differences in CDFs{p_end}

{dlgtab:Graph Appearance}

{phang}
{opt title(string)} title for the graph. Defaults vary by plot type.

{phang}
{opt ytitle(string)} y-axis title. Defaults vary by plot type.

{phang}
{opt xtitle(string)} x-axis title. Defaults vary by plot type.

{phang}
{opt color1(string)} color for first series (default: "blue").

{phang}
{opt color2(string)} color for second series in level plots (default: "red").

{phang}
{opt cicolor(string)} color for confidence intervals (default: "gs12").

{phang}
{opt lwidth(string)} line width (default: "medium").

{phang}
{opt lpattern(string)} line pattern for second series (default: "dash").

{phang}
{opt legend(string)} legend options. See {help legend_options}.

{phang}
{opt byopts(string)} options for by() graphs. See {help by_option}.

{phang}
{opt plotregion(string)} plot region options. See {help region_options}.

{phang}
{opt graphregion(string)} graph region options. See {help region_options}.

{phang}
{opt scheme(string)} scheme name. See {help scheme_option}.

{marker examples}
{title:Examples}

{pstd}Basic quantile difference plot after disco estimation:{p_end}
{phang2}{cmd:. disco y id time, idtarget(1) t0(10)}{p_end}
{phang2}{cmd:. disco_plot}{p_end}

{pstd}CDF plot with custom styling:{p_end}
{phang2}{cmd:. disco y id time, idtarget(1) t0(10) agg("cdf")}{p_end}
{phang2}{cmd:. disco_plot, color1(navy) color2(maroon) lwidth(thick)}{p_end}

{pstd}Quantile plot with confidence intervals and custom title:{p_end}
{phang2}{cmd:. disco y id time, idtarget(1) t0(10) ci}{p_end}
{phang2}{cmd:. disco_plot, agg(quantile) title("Distribution Effects Over Time")}{p_end}

{pstd}CDF differences with custom legend:{p_end}
{phang2}{cmd:. disco y id time, idtarget(1) t0(10) agg("cdfDiff")}{p_end}
{phang2}{cmd:. disco_plot, legend(ring(0) pos(11) rows(2))}{p_end}

{marker details}
{title:Details}

{pstd}
The command automatically handles different types of plots based on the aggregation type specified 
either in the original {cmd:disco} command or through the {cmd:agg()} option in {cmd:disco_plot}. 
The default visualization changes based on this type:

{phang2}For {cmd:quantileDiff} and {cmd:cdfDiff}:{p_end}
{phang3}- Shows differences between treated and synthetic units{p_end}
{phang3}- Includes confidence intervals if CI was specified in {cmd:disco}{p_end}
{phang3}- Uses single line plots with optional confidence bands{p_end}

{phang2}For {cmd:quantile} and {cmd:cdf}:{p_end}
{phang3}- Shows levels for both treated and synthetic units{p_end}
{phang3}- Uses dual line plots with different patterns{p_end}
{phang3}- Includes appropriate legends{p_end}

{pstd}
All plots are created as small multiples using Stata's {cmd:by()} functionality, with one panel 
per time period. This allows for easy visualization of how distributional effects evolve over time.

{marker results}
{title:Stored results}

{pstd}
{cmd:disco_plot} does not store results but creates graphs using the results stored by {cmd:disco}.

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