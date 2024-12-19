program define disco_plot
    version 18.0
    
    syntax, AGG(string) M(integer) G(integer) T_max(integer) DOCI(integer) CL(real)
    
    preserve
    clear
    
    local cl_txt = subinstr("`cl'", ".", "", .)

    if "`agg'" == "quantileDiff" {
        local ytitle "Difference in Quantile Functions"
		local title "Distributional Effects by Time Period"


        if `doci' {
            set obs `m'
            gen tau = (_n-1)/(`m'-1)
            
            svmat quantile_diff
            svmat qdiff_lower
            svmat qdiff_upper
            reshape long quantile_diff qdiff_lower qdiff_upper, i(tau) j(time)
            
            twoway (rarea qdiff_lower qdiff_upper tau, color(gs12)) ///
                   (line quantile_diff tau, lcolor(black)), ///
                   by(time, note("") title("`title'")) ///
                   ytitle("`ytitle'") ///
                   xtitle("Quantile") legend(off) ///
                   legend(label(1 "`cl_txt' % confidence intervals") label(2 "Estimates")) ///
                   xlabel(0(.2)1) ylabel(, angle(horizontal))
        }
        else {
            set obs `m'
            gen tau = (_n-1)/(`m'-1)
            
            svmat quantile_diff
            reshape long quantile_diff, i(tau) j(time)
            
            twoway line quantile_diff tau, ///
                   by(time, note("")) ///
                   title("`title'") ytitle("`ytitle'") ///
                   xtitle("Quantile") legend(off) ///
                   xlabel(0(.2)1) ylabel(, angle(horizontal))
        }
    }
    else if "`agg'" == "quantile" {
        local ytitle "Quantile Function (Synthetic vs. Target)"
        local title "Synthetic vs. Treated Quantiles by Time Period"
        
        set obs `m'
        gen tau = (_n-1)/(`m'-1)
        
        svmat quantile_t
        svmat quantile_synth
        reshape long quantile_t quantile_synth, i(tau) j(time)
        

        twoway (line quantile_t tau, lcolor(blue)) ///
               (line quantile_synth tau, lcolor(red) lpattern(dash)), ///
               by(time, note("") title("`title'") ) ///
               ytitle("`ytitle'") ///
               xtitle("Quantile") ///
               legend(order(1 "Observed" 2 "Synthetic") ring(0) pos(1)) ///
               xlabel(0(.2)1) ylabel(, angle(horizontal))
    }
    else if inlist("`agg'", "cdf", "cdfDiff") {
        if "`agg'" == "cdfDiff" {
            local ytitle "Difference in CDFs"
        }
        else {
            local ytitle "CDF (Synthetic vs. Target)"
        }
        
        local gmin = e(amin)
        local gmax = e(amax)
        
        if "`agg'" == "cdfDiff" {
			local title "Distributional Effects by Time Period"

            if `doci' {
                set obs `g'
                gen grid_val = `gmin' + (_n-1)*(`gmax' - `gmin')/(`g'-1)
        
                svmat cdf_diff 
                svmat cdiff_lower
                svmat cdiff_upper
                reshape long cdf_diff cdiff_lower cdiff_upper, i(grid_val) j(time)
                

                
                twoway (rarea cdiff_lower cdiff_upper grid_val, color(gs12)) ///
                       (line cdf_diff grid_val, lcolor(black)), ///
                       by(time, note("") title("`title'")) ///
                       legend(label(1 "`cl_txt' % confidence intervals") label(2 "Estimates")) ///
                       ytitle("`ytitle'") ///
                       xtitle("Y") legend(off) ///
                       ylabel(, angle(horizontal))
            }
            else {
                set obs `g'
                gen grid_val = `gmin' + (_n-1)*(`gmax' - `gmin')/(`g'-1)
                
                svmat cdf_diff
                reshape long cdf_diff, i(grid_val) j(time)

            
                twoway line cdf_diff grid_val, ///
                       by(time, note("") title("`title'")) ///
                       ytitle("`ytitle'") ///
                       xtitle("Y") legend(off) ///
                       ylabel(, angle(horizontal))
            }
        }
        else {
			local title = "Synthetic vs. Treated CDFs by Time Period"
            set obs `g'
            gen grid_val = `gmin' + (_n-1)*(`gmax' - `gmin')/(`g'-1)
            svmat cdf_t
            svmat cdf_synth
            reshape long cdf_t cdf_synth, i(grid_val) j(time)
            
            twoway (line cdf_t grid_val, lcolor(blue)) ///
                   (line cdf_synth grid_val, lcolor(red) lpattern(dash)), ///
                   by(time, note("") title("`title'")) ///
                   ytitle("`ytitle'") ///
                   xtitle("Y") ///
                   legend(order(1 "Observed" 2 "Synthetic") ring(0) pos(1)) ///
                   ylabel(0(.2)1, angle(horizontal))
        }
    }
    
    restore
end
