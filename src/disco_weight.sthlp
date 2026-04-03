{smcl}
{* *! version 1.0.0 19dec2024}{...}
{vieweralsosee "[ST] disco" "help disco"}{...}
{vieweralsosee "[ST] disco_estat" "help disco_estat"}{...}
{vieweralsosee "[ST] disco_plot" "help disco_plot"}{...}
{viewerjumpto "Syntax" "disco_weight##syntax"}{...}
{viewerjumpto "Description" "disco_weight##description"}{...}
{viewerjumpto "Options" "disco_weight##options"}{...}
{viewerjumpto "Examples" "disco_weight##examples"}{...}
{viewerjumpto "Stored results" "disco_weight##results"}{...}

{title:disco_weight}

{phang}
{bf:disco_weight} {hline 2} Post-estimation command to display and store synthetic control weights with unit names

{marker description}
{title:Description}

{pstd}
{cmd:disco_weight} displays and stores the top synthetic control weights after {cmd:disco} estimation, 
matching the numeric unit IDs with their corresponding names. Results are stored in a new frame 
for further analysis.

{marker syntax}
{title:Syntax}

{p 8 17 2}
{cmdab:disco_weight} {it:id_var} {it:name_var} [{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Options}
{synopt:{opt n(#)}}number of top weights to display (default: 5){p_end}
{synopt:{opt f:ormat(str)}}display format for weights (default: %12.4f){p_end}
{synopt:{opt frame(str)}}name of frame to store results (default: disco_weights){p_end}
{synopt:{opt r:ound(#)}}rounding precision for weights (default: 0.0001){p_end}
{synoptline}

{pstd}
{cmd:disco_weight} requires that {cmd:disco} has been run previously.

{marker options}
{title:Options}

{dlgtab:Main}

{phang}
{opt n(#)} specifies the number of top weights to display and store. The default is 5.

{phang}
{opt format(str)} specifies the display format for weights. The default is %12.4f.

{phang}
{opt frame(str)} specifies the name of the frame where results will be stored. 
The default is "disco_weights". If a frame with this name already exists, it will be replaced.

{phang}
{opt round(#)} specifies the rounding precision for weights. The default is 0.0001.

{marker examples}
{title:Examples}

{pstd}Display top 5 weights with company names after running disco:{p_end}
{phang2}{cmd:. disco_weight id_col company_name}{p_end}

{pstd}Store top 10 weights in a custom frame:{p_end}
{phang2}{cmd:. disco_weight id_col company_name, n(10) frame(my_weights)}{p_end}

{pstd}Access the stored weights:{p_end}
{phang2}{cmd:. frame my_weights: list name weight}{p_end}

{marker results}
{title:Stored results}

{pstd}
{cmd:disco_weight} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n)}}number of top weights displayed{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}"disco_weight"{p_end}
{synopt:{cmd:r(format)}}display format used{p_end}
{synopt:{cmd:r(frame)}}name of frame where results are stored{p_end}

{pstd}
Results are also stored in a new frame (specified by {cmd:frame()} option) with the following variables:

{p2col 5 20 24 2: Variables}{p_end}
{synopt:{cmd:name}}unit names corresponding to weights{p_end}
{synopt:{cmd:weight}}synthetic control weights{p_end}

{marker author}
{title:Author}

{pstd}
David Van Dijcke{break}
University of Michigan, Ann Arbor{break}
{browse "mailto:dvdijcke@umich.edu":dvdijcke@umich.edu}

{title:Version}

{pstd}
1.0.0 (December 2024)
{p_end}