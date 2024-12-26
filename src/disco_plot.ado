
program define disco_plot
    version 18.0
    
    syntax, AGG(string) M(integer) G(integer) T_max(integer) DOCI(integer) CL(real) ///
        quantile_diff(name) quantile_t(name) quantile_synth(name) ///
        cdf_diff(name) cdf_synth(name) cdf_t(name) ///
        [qdiff_lower(name) qdiff_upper(name) cdiff_lower(name) cdiff_upper(name) ///
        TITLE(string) YTITLE(string) XTITLE(string) ///
        BYOPTS(string) PLOTREGION(string) GRAPHREGION(string) ///
        SCHEME(string) LEGEND(string) COLOR1(string) COLOR2(string) ///
        CIcolor(string) LWIDTH(string) LPATTERN(string)]
    
    preserve
    clear
    
    // Set default colors and styles if not specified
    if "`color1'" == "" local color1 "blue"
    if "`color2'" == "" local color2 "red" 
    if "`cicolor'" == "" local cicolor "gs12"
    if "`lwidth'" == "" local lwidth "medium"
    if "`lpattern'" == "" local lpattern "dash"
    
    local cl_txt = subinstr("`cl'", ".", "", .)

    if "`agg'" == "quantileDiff" {
        // Set default titles if not specified
        if "`title'" == "" local title "Distributional Effects by Time Period"
        if "`ytitle'" == "" local ytitle "Difference in Quantile Functions"
        if "`xtitle'" == "" local xtitle "Quantile"

        if `doci' {
            set obs `m'
            gen tau = (_n-1)/(`m'-1)
            
            svmat `quantile_diff'
            svmat `qdiff_lower'
            svmat `qdiff_upper'
            reshape long `quantile_diff' `qdiff_lower' `qdiff_upper', i(tau) j(time)
            
            twoway (rarea `qdiff_lower' `qdiff_upper' tau, color(`cicolor')) ///
                   (line `quantile_diff' tau, lcolor(`color1') lwidth(`lwidth')), ///
                   by(time, note("") title(`"`title'"') `byopts') ///
                   ytitle(`"`ytitle'"') xtitle(`"`xtitle'"') ///
                   xlabel(0(.2)1) ylabel(, angle(horizontal)) ///
                   legend(label(1 "`cl_txt'% confidence intervals") ///
                          label(2 "Estimates") `legend') ///
                   `plotregion' `graphregion' `scheme'
        }
        else {
            set obs `m'
            gen tau = (_n-1)/(`m'-1)
            
            svmat `quantile_diff'
            reshape long `quantile_diff', i(tau) j(time)
            
            twoway line `quantile_diff' tau, ///
                   lcolor(`color1') lwidth(`lwidth') ///
                   by(time, note("") title(`"`title'"') `byopts') ///
                   ytitle(`"`ytitle'"') xtitle(`"`xtitle'"') ///
                   xlabel(0(.2)1) ylabel(, angle(horizontal)) ///
                   legend(off) ///
                   `plotregion' `graphregion' `scheme'
        }
    }
    else if "`agg'" == "quantile" {
        // Set default titles
        if "`title'" == "" local title "Synthetic vs. Treated Quantiles by Time Period"
        if "`ytitle'" == "" local ytitle "Quantile Function (Synthetic vs. Target)"
        if "`xtitle'" == "" local xtitle "Quantile"
        
        set obs `m'
        gen tau = (_n-1)/(`m'-1)
        
        svmat `quantile_t'
        svmat `quantile_synth'
        reshape long `quantile_t' `quantile_synth', i(tau) j(time)
        
        twoway (line `quantile_t' tau, lcolor(`color1') lwidth(`lwidth')) ///
               (line `quantile_synth' tau, lcolor(`color2') lwidth(`lwidth') lpattern(`lpattern')), ///
               by(time, note("") title(`"`title'"') `byopts') ///
               ytitle(`"`ytitle'"') xtitle(`"`xtitle'"') ///
               legend(order(1 "Observed" 2 "Synthetic") ring(0) pos(1) `legend') ///
               xlabel(0(.2)1) ylabel(, angle(horizontal)) ///
               `plotregion' `graphregion' `scheme'
    }
    else if inlist("`agg'", "cdf", "cdfDiff") {
        // Set default titles
        if "`agg'" == "cdfDiff" {
            if "`ytitle'" == "" local ytitle "Difference in CDFs"
            if "`title'" == "" local title "Distributional Effects by Time Period"
        }
        else {
            if "`ytitle'" == "" local ytitle "CDF (Synthetic vs. Target)"
            if "`title'" == "" local title "Synthetic vs. Treated CDFs by Time Period"
        }
        if "`xtitle'" == "" local xtitle "Y"
        
        local gmin = `min'  // You'll need to pass these from disco.ado TODO: this is a bug
        local gmax = `max'
        
        if "`agg'" == "cdfDiff" {
            if `doci' {
                set obs `g'
                gen grid_val = `gmin' + (_n-1)*(`gmax' - `gmin')/(`g'-1)
        
                svmat `cdf_diff'
                svmat `cdiff_lower'
                svmat `cdiff_upper'
                reshape long `cdf_diff' `cdiff_lower' `cdiff_upper', i(grid_val) j(time)
                
                twoway (rarea `cdiff_lower' `cdiff_upper' grid_val, color(`cicolor')) ///
                       (line `cdf_diff' grid_val, lcolor(`color1') lwidth(`lwidth')), ///
                       by(time, note("") title(`"`title'"') `byopts') ///
                       ytitle(`"`ytitle'"') xtitle(`"`xtitle'"') ///
                       legend(label(1 "`cl_txt'% confidence intervals") ///
                              label(2 "Estimates") `legend') ///
                       ylabel(, angle(horizontal)) ///
                       `plotregion' `graphregion' `scheme'
            }
            else {
                set obs `g'
                gen grid_val = `gmin' + (_n-1)*(`gmax' - `gmin')/(`g'-1)
                
                svmat `cdf_diff'
                reshape long `cdf_diff', i(grid_val) j(time)
                
                twoway line `cdf_diff' grid_val, ///
                       lcolor(`color1') lwidth(`lwidth') ///
                       by(time, note("") title(`"`title'"') `byopts') ///
                       ytitle(`"`ytitle'"') xtitle(`"`xtitle'"') ///
                       legend(off) ylabel(, angle(horizontal)) ///
                       `plotregion' `graphregion' `scheme'
            }
        }
        else {
            set obs `g'
            gen grid_val = `gmin' + (_n-1)*(`gmax' - `gmin')/(`g'-1)
            
            svmat `cdf_t'
            svmat `cdf_synth'
            reshape long `cdf_t' `cdf_synth', i(grid_val) j(time)
            
            twoway (line `cdf_t' grid_val, lcolor(`color1') lwidth(`lwidth')) ///
                   (line `cdf_synth' grid_val, lcolor(`color2') lwidth(`lwidth') lpattern(`lpattern')), ///
                   by(time, note("") title(`"`title'"') `byopts') ///
                   ytitle(`"`ytitle'"') xtitle(`"`xtitle'"') ///
                   legend(order(1 "Observed" 2 "Synthetic") ring(0) pos(1) `legend') ///
                   ylabel(0(.2)1, angle(horizontal)) ///
                   `plotregion' `graphregion' `scheme'
        }
    }
    
    restore
end
