

// Modified disco_estat.ado
program define disco_estat, rclass
    version 16.0
    
    if "`e(cmd)'" != "disco" {
        error 301
    }
    
    if "`0'" == "summary" {
        if "`e(agg)'" == "" {
            di as error "No aggregation statistics were computed. Rerun disco with agg() option."
            exit 198
        }
        
        if "`e(agg)'" == "quantileDiff" {
            local agtxt = "quantile"
            local fmt "%4.2f"  // Format for probabilities [0,1]
        } 
        else if "`e(agg)'" == "cdfDiff" {
            local agtxt = "CDF"
            local fmt "%9.2f"  // Format for Y values
        }
        
        tempname stats
        matrix `stats' = e(summary_stats)
        
        // Get number of rows
        local nr = rowsof(`stats')
        
        // Display header
        di _n as txt "Summary of `agtxt' effects"
        di as txt "{hline 80}"
        di as txt "Time period   Range                    Effect     Std. Err.    [`e(cl)'% Conf. Interval]"
        di as txt "{hline 80}"
        
        // Display results
        forvalues i = 1/`nr' {
            local t = `stats'[`i',1]
            local qstart = `stats'[`i',2]
            local qend = `stats'[`i',3]
            local effect = `stats'[`i',4]
            local se = `stats'[`i',5]
            local ci_l = `stats'[`i',6]
            local ci_u = `stats'[`i',7]
            
            // Check if significant (CI doesn't include 0)
            local sig = (`ci_l' > 0 | `ci_u' < 0) * "*"
            
            di as txt %9.0g `t' "    " ///
               as txt `fmt' `qstart' "-" `fmt' `qend' "    " ///
               as res %9.3f `effect' "    " ///
               as res %9.3f `se' "    " ///
               as res %9.3f `ci_l' "  " %9.3f `ci_u' ///
               as txt "`sig'"
        }
        
        di as txt "{hline 80}"
        di as txt "* denotes significance at `e(cl)'% confidence level"
    }
    else {
        di as error "unknown subcommand"
        exit 198
    }
end
