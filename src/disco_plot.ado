/*******************************************************************************
 * Distributional Synthetic Controls Plot (disco_plot)
 * 
 * Post-estimation command for DiSCo (disco.ado).
 * Reads results from e(...) and plots them by time period.
 *
 * Implementation based on Gunsilius (2023) and Van Dijcke, Gunsilius, and
 * Wright (2024).
 *
 * Syntax (after disco):
 *     disco_plot [ , AGG(string) M(integer) G(integer) DOCI(integer) 
 *                      CL(real) TITLE(string) YTITLE(string) XTITLE(string) 
 *                      COLOR1(string) COLOR2(string) CIcolor(string) 
 *                      LWIDTH(string) LPATTERN(string) BYOPTS(string) 
 *                      PLOTREGION(string) GRAPHREGION(string) SCHEME(string) 
 *                      LEGEND(string)
 *                ]
 *
 * Author: David Van Dijcke
 * Version: 1.0.0
 * Date: December 2024
 *******************************************************************************/
program define disco_plot, eclass
    version 18.0

    // 2. Syntax parsing
    syntax[ , AGG(string) /// 
            TITLE(string) YTITLE(string) XTITLE(string) ///
            COLOR1(string) COLOR2(string) CIcolor(string) ///
            LWIDTH(string) LPATTERN(string) LEGEND(string) ///
            BYOPTS(string) PLOTREGION(string) GRAPHREGION(string) ///
            SCHEME(string)]
		
	// 1. Check that e(cmd) is "disco"
    if "`e(cmd)'" != "disco" {
        di as error "disco_plot is a post-estimation command for disco. Run disco first."
        exit 198
    }
    // set default agg if user didnt specify
    if missing(e(agg)) {
		local agg = "quantileDiff"
	} 
	else {
		local agg = e(agg)
	}
    local m = e(m)
    local g = e(g)
	local doci = e(doci)
	local cl = e(cl)
	

    // Confidence intervals (if e(doci)==1 by default)
    if missing(`doci') local doci = e(doci)
    if missing(`cl') local cl = e(cl)
    local cl_txt = subinstr("`cl'", ".", "", .)

    // Extract needed matrices from e()
    matrix quantile_diff    = e(quantile_diff)
    matrix quantile_t       = e(quantile_t)
    matrix quantile_synth   = e(quantile_synth)
    matrix cdf_diff         = e(cdf_diff)
    matrix cdf_synth        = e(cdf_synth)
    matrix cdf_t            = e(cdf_t)

    // If we have CIs
    if `doci' == 1 {
        matrix qdiff_lower   = e(qdiff_lower)
        matrix qdiff_upper   = e(qdiff_upper)
        matrix cdiff_lower   = e(cdiff_lower)
        matrix cdiff_upper   = e(cdiff_upper)
    }

    // Scalars for time range and x-lims
    scalar t_max = e(t_max)
    scalar xmin  = e(amin)
    scalar xmax  = e(amax)

    // Default graph options
    if "`color1'" == "" local color1 "blue"
    if "`color2'" == "" local color2 "red"
    if "`cicolor'" == "" local cicolor "gs12"
    if "`lwidth'" == "" local lwidth "medium"
    if ("`lpattern'" == "") local lpattern "dash"

    // Now replicate your plotting logic using only the e(...) results.
    preserve
    clear
    
    if "`agg'" == "quantileDiff" {
        // Default titles
        if "`title'" == "" local title "Distributional Effects by Time Period"
        if "`ytitle'" == "" local ytitle "Difference in Quantile Functions"
        if "`xtitle'" == "" local xtitle "Quantile"

        if `doci' {
            quietly: set obs `m'
            gen tau = (_n-1)/(`m'-1)
            
            svmat quantile_diff
            svmat qdiff_lower
            svmat qdiff_upper
            quietly: reshape long quantile_diff qdiff_lower qdiff_upper, i(tau) j(time)
            
            twoway (rarea qdiff_lower qdiff_upper tau, color(`cicolor')) ///
                   (line quantile_diff tau, lcolor(`color1') lwidth(`lwidth')), ///
                   by(time, note("") title("`title'") `byopts') ///
                   ytitle("`ytitle'") xtitle("`xtitle'") ///
                   xlabel(0(.2)1) ylabel(, angle(horizontal)) ///
                   legend(label(1 "`cl_txt'% CIs") label(2 "Estimates") `legend') ///
                   `plotregion' `graphregion' `scheme'
        }
        else {
            quietly: set obs `m'
            gen tau = (_n-1)/(`m'-1)
            
            svmat quantile_diff
            quietly: reshape long quantile_diff, i(tau) j(time)
            
            twoway line quantile_diff tau, ///
                   lcolor(`color1') lwidth(`lwidth') ///
                   by(time, note("") title("`title'") `byopts') ///
                   ytitle("`ytitle'") xtitle("`xtitle'") ///
                   xlabel(0(.2)1) ylabel(, angle(horizontal)) ///
                   legend(off) ///
                   `plotregion' `graphregion' `scheme'
        }
    }
    else if "`agg'" == "quantile" {
        // Default titles
        if "`title'" == "" local title "Synthetic vs. Treated Quantiles by Time Period"
        if "`ytitle'" == "" local ytitle "Quantile Function (Synthetic vs. Target)"
        if "`xtitle'" == "" local xtitle "Quantile"
        
        quietly: set obs `m'
        gen tau = (_n-1)/(`m'-1)
        
        svmat quantile_t
        svmat quantile_synth
        quietly: reshape long quantile_t quantile_synth, i(tau) j(time)
        
        twoway (line quantile_t tau, lcolor(`color1') lwidth(`lwidth')) ///
               (line quantile_synth tau, lcolor(`color2') lwidth(`lwidth') lpattern(`lpattern')), ///
               by(time, note("") title("`title'") `byopts') ///
               ytitle("`ytitle'") xtitle("`xtitle'") ///
               legend(order(1 "Observed" 2 "Synthetic") ring(0) pos(1) `legend') ///
               xlabel(0(.2)1) ylabel(, angle(horizontal)) ///
               `plotregion' `graphregion' `scheme'
    }
    else if inlist("`agg'", "cdf", "cdfDiff") {
        // Defaults
        if "`agg'" == "cdfDiff" {
            if "`ytitle'" == "" local ytitle "Difference in CDFs"
            if "`title'" == "" local title "Distributional Effects by Time Period"
        }
        else {
            if "`ytitle'" == "" local ytitle "CDF (Synthetic vs. Target)"
            if "`title'" == "" local title "Synthetic vs. Treated CDFs by Time Period"
        }
        if "`xtitle'" == "" local xtitle "Y"
        
        // cdfDiff
        if "`agg'" == "cdfDiff" {
            if `doci' {
                quietly: set obs `g'
                quietly: gen grid_val = xmin + (_n-1)*(xmax - xmin)/(`g'-1)
        
                svmat cdf_diff
                svmat cdiff_lower
                svmat cdiff_upper
                quietly: reshape long cdf_diff cdiff_lower cdiff_upper, i(grid_val) j(time)
                
                twoway (rarea cdiff_lower cdiff_upper grid_val, color(`cicolor')) ///
                       (line cdf_diff grid_val, lcolor(`color1') lwidth(`lwidth')), ///
                       by(time, note("") title("`title'") `byopts') ///
                       ytitle("`ytitle'") xtitle("`xtitle'") ///
                       legend(label(1 "`cl_txt'% CIs") label(2 "Estimates") `legend') ///
                       ylabel(, angle(horizontal)) ///
                       `plotregion' `graphregion' `scheme'
            }
            else {
                quietly: set obs `g'
                gen grid_val = xmin + (_n-1)*(xmax - xmin)/(`g'-1)
                
                svmat cdf_diff
                quietly: reshape long cdf_diff, i(grid_val) j(time)
                
                twoway line cdf_diff grid_val, ///
                       lcolor(`color1') lwidth(`lwidth') ///
                       by(time, note("") title("`title'") `byopts') ///
                       ytitle("`ytitle'") xtitle("`xtitle'") ///
                       legend(off) ylabel(, angle(horizontal)) ///
                       `plotregion' `graphregion' `scheme'
            }
        }
        // cdf (levels)
        else {
            quietly: set obs `g'
            gen grid_val = xmin + (_n-1)*(xmax - xmin)/(`g'-1)
            
            svmat cdf_t
            svmat cdf_synth
            quietly: reshape long cdf_t cdf_synth, i(grid_val) j(time)
            
            twoway (line cdf_t grid_val, lcolor(`color1') lwidth(`lwidth')) ///
                   (line cdf_synth grid_val, lcolor(`color2') lwidth(`lwidth') lpattern(`lpattern')), ///
                   by(time, note("") title("`title'") `byopts') ///
                   ytitle("`ytitle'") xtitle("`xtitle'") ///
                   legend(order(1 "Observed" 2 "Synthetic") ring(0) pos(1) `legend') ///
                   ylabel(0(.2)1, angle(horizontal)) ///
                   `plotregion' `graphregion' `scheme'
        }
    }

    restore
end
