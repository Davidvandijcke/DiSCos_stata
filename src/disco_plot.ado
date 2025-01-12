program define disco_plot, eclass
    version 18.0

    // 2. Syntax parsing
    syntax[ , AGG(string) CATEGORICAL /// Add categorical option
            TITLE(string) YTITLE(string) XTITLE(string) ///
            COLOR1(string) COLOR2(string) CIcolor(string) ///
            LWIDTH(string) LPATTERN(string) LEGEND(string) ///
            BYOPTS(string) PLOTREGION(string) GRAPHREGION(string) ///
            SCHEME(string) HLINE(real 999) VLINE(real 999) ///
            XRANGE(numlist min=2 max=2) YRANGE(numlist min=2 max=2)]
		
    // 1. Check that e(cmd) is "disco"
    if "`e(cmd)'" != "disco" {
        di as error "disco_plot is a post-estimation command for disco. Run disco first."
        exit 198
    }
    
    // set default agg if user didnt specify
    if "`agg'" == "" {
        if missing(e(agg)) {
            local agg = "quantileDiff"
        } 
        else {
            local agg = e(agg)
        }
    }
    
    local g = e(g)
    local doci = e(doci)
    local cl = e(cl)

    // Initialize line options as empty
    local hline_opt ""
    local vline_opt ""
    local xrange_opt ""
    local yrange_opt ""
    
    // Add horizontal line if specified
    if `hline' != 999 {
        local hline_opt = "yline(`hline', lcolor(gray) lpattern(dash))"
    }
    
    // Add vertical line if specified
    if `vline' != 999 {
        local vline_opt "xline(`vline', lcolor(gray) lpattern(dash))"
    }
    
    // Handle axis ranges if specified
    if "`xrange'" != "" {
        tokenize `xrange'
        local xrange_opt "xscale(range(`1' `2'))"
    }
    
    if "`yrange'" != "" {
        tokenize `yrange'
        local yrange_opt "yscale(range(`1' `2'))"
    }

    // Confidence intervals (if e(doci)==1 by default)
    if missing(`doci') local doci = e(doci)
    if missing(`cl') local cl = e(cl)
    // Convert confidence level from decimal to percentage (e.g., 0.9 -> 90)
    local cl_num = real("`cl'")*100
    local cl_txt = string(`cl_num', "%3.0f")

    // Extract needed matrices from e()
    matrix quantile_diff = e(quantile_diff)
    matrix quantile_t = e(quantile_t)
    matrix quantile_synth = e(quantile_synth)
    matrix cdf_diff = e(cdf_diff)
    matrix cdf_synth = e(cdf_synth)
    matrix cdf_t = e(cdf_t)

    // If we have CIs
    if `doci' == 1 {
        matrix qdiff_lower = e(qdiff_lower)
        matrix qdiff_upper = e(qdiff_upper)
        matrix cdiff_lower = e(cdiff_lower)
        matrix cdiff_upper = e(cdiff_upper)
    }

    // Scalars for time range and x-lims
    scalar t_max = e(t_max)
    scalar xmin = e(amin)
    scalar xmax = e(amax)

    // Default graph options
    if "`color1'" == "" local color1 "blue"
    if "`color2'" == "" local color2 "red"
    if "`cicolor'" == "" local cicolor "gs12"
    if "`lwidth'" == "" local lwidth "medium"
    if ("`lpattern'" == "") local lpattern "dash"

    preserve
    clear
    
    // Handle standard quantileDiff and quantile cases as before
    if "`agg'" == "quantileDiff" {
        // Default titles
        if "`title'" == "" local title "Distributional Effects by Time Period"
        if "`ytitle'" == "" local ytitle "Difference in Quantile Functions"
        if "`xtitle'" == "" local xtitle "Quantile"

        if `doci' {
            quietly: set obs `g'
            gen tau = (_n-1)/(`g'-1)
            
            svmat quantile_diff
            svmat qdiff_lower
            svmat qdiff_upper
            quietly: reshape long quantile_diff qdiff_lower qdiff_upper, i(tau) j(time)
            
            twoway (rarea qdiff_lower qdiff_upper tau, color(`cicolor')) ///
                   (line quantile_diff tau, lcolor(`color1') lwidth(`lwidth')), ///
                   `hline_opt' `vline_opt' ///
                   by(time, note("") title("`title'") `byopts') ///
                   ytitle("`ytitle'") xtitle("`xtitle'") ///
                   `xrange_opt' `yrange_opt' ///
                   legend(label(1 "`cl_txt'% CIs") label(2 "Estimates") `legend') ///
                   plotregion(`plotregion') graphregion(`graphregion') scheme(`scheme')
        }
        else {
            quietly: set obs `g'
            gen tau = (_n-1)/(`g'-1)
            
            svmat quantile_diff
            quietly: reshape long quantile_diff, i(tau) j(time)
            
            twoway line quantile_diff tau, ///
                   `hline_opt' `vline_opt' ///
                   lcolor(`color1') lwidth(`lwidth') ///
                   by(time, note("") title("`title'") `byopts') ///
                   ytitle("`ytitle'") xtitle("`xtitle'") ///
                   `xrange_opt' `yrange_opt' ///
                   legend(off) ///
                   plotregion(`plotregion') graphregion(`graphregion') scheme(`scheme')
        }
    }
    else if "`agg'" == "quantile" {
        // Default titles
        if "`title'" == "" local title "Synthetic vs. Treated Quantiles by Time Period"
        if "`ytitle'" == "" local ytitle "Quantile Function (Synthetic vs. Target)"
        if "`xtitle'" == "" local xtitle "Quantile"
        
        quietly: set obs `g'
        gen tau = (_n-1)/(`g'-1)
        
        svmat quantile_t
        svmat quantile_synth
        quietly: reshape long quantile_t quantile_synth, i(tau) j(time)
        
        twoway (line quantile_t tau, lcolor(`color1') lwidth(`lwidth')) ///
               (line quantile_synth tau, lcolor(`color2') lwidth(`lwidth') lpattern(`lpattern')), ///
               `hline_opt' `vline_opt' ///
               by(time, note("") title("`title'") `byopts') ///
               ytitle("`ytitle'") xtitle("`xtitle'") ///
               `xrange_opt' `yrange_opt' ///
               legend(order(1 "Observed" 2 "Synthetic") ring(0) pos(1) `legend') ///
               plotregion(`plotregion') graphregion(`graphregion') scheme(`scheme')
    }
    // Modified CDF section to handle categorical option
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
        
        // Categorical plot for cdfDiff
        if "`agg'" == "cdfDiff" & "`categorical'" != "" {
            if `doci' {
                quietly: set obs `g'
                quietly: gen category = _n
                
                svmat cdf_diff
                svmat cdiff_lower
                svmat cdiff_upper
                quietly: reshape long cdf_diff cdiff_lower cdiff_upper, i(category) j(time)
                
                // Create x-axis labels based on the range in the data
                quiet summ category
                local min = r(min)
                local max = r(max)
                local xlabel_opt "xlabel(`min'(1)`max', angle(45))"
                
                // Create categorical bar plot with CIs
                twoway (bar cdf_diff category, color(`color1') barwidth(0.8)) ///
                       (rcap cdiff_upper cdiff_lower category, color(navy)), ///
                       `hline_opt' `vline_opt' ///
                       `xrange_opt' `yrange_opt' ///
                       by(time, note("") title("`title'") `byopts') ///
                       ytitle("`ytitle'") xtitle("`xtitle'") ///
                       `xlabel_opt' ///
                       legend(label(1 "Effect") label(2 "`cl_txt'% CIs") `legend') ///
                       plotregion(`plotregion') graphregion(`graphregion') scheme(`scheme')
            }
            else {
                quietly: set obs `g'
                gen category = _n
                
                svmat cdf_diff
                quietly: reshape long cdf_diff, i(category) j(time)
                
                // Get category labels from matrix column names if available
                local colnames : colnames cdf_diff
                if "`colnames'" != "" {
                    gen category_label = ""
                    local i = 1
                    foreach name of local colnames {
                        quietly replace category_label = "`name'" if category == `i'
                        local i = `i' + 1
                    }
                }
                else {
                    // If no column names, use numeric labels
                    gen category_label = string(category)
                }
                
                // Create categorical bar plot without CIs
                twoway bar cdf_diff category, ///
                       `hline_opt' `vline_opt' ///
                       `xrange_opt' `yrange_opt' ///
                       color(`color1') barwidth(0.8) ///
                       by(time, note("") title("`title'") `byopts') ///
                       ytitle("`ytitle'") xtitle("`xtitle'") ///
                       xlabel(1(1)10, valuelabel angle(45)) ///
                       legend(off) ///
                       plotregion(`plotregion') graphregion(`graphregion') scheme(`scheme')
            }
        }
        // Regular CDF plots (non-categorical)
        else {
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
                           `hline_opt' `vline_opt' ///
                           `xrange_opt' `yrange_opt' ///
                           by(time, note("") title("`title'") `byopts') ///
                           ytitle("`ytitle'") xtitle("`xtitle'") ///
                           legend(label(1 "`cl_txt'% CIs") label(2 "Estimates") `legend') ///
                           plotregion(`plotregion') graphregion(`graphregion') scheme(`scheme')
                }
                else {
                    quietly: set obs `g'
                    gen grid_val = xmin + (_n-1)*(xmax - xmin)/(`g'-1)
                    
                    svmat cdf_diff
                    quietly: reshape long cdf_diff, i(grid_val) j(time)
                    
                    twoway line cdf_diff grid_val, ///
                           `hline_opt' `vline_opt' ///
                           `xrange_opt' `yrange_opt' ///
                           lcolor(`color1') lwidth(`lwidth') ///
                           by(time, note("") title("`title'") `byopts') ///
                           ytitle("`ytitle'") xtitle("`xtitle'") ///
                           legend(off) ///
                           plotregion(`plotregion') graphregion(`graphregion') scheme(`scheme')
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
                       `hline_opt' `vline_opt' ///
                       `xrange_opt' `yrange_opt' ///
                       by(time, note("") title("`title'") `byopts') ///
                       ytitle("`ytitle'") xtitle("`xtitle'") ///
                       legend(order(1 "Observed" 2 "Synthetic") ring(0) pos(1) `legend') ///
                       plotregion(`plotregion') graphregion(`graphregion') scheme(`scheme')
            }
        }
    }

    restore
end
