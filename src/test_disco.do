clear
set seed 123

capture program drop _all
mata: mata clear

net install disco, from("https://raw.githubusercontent.com/Davidvandijcke/DiSCos_stata/main/src/") replace
//net install disco, from("/Users/davidvandijcke/University of Michigan Dropbox/David Van Dijcke/Flo_GSRA/stata_repo/src/") replace
// do disco.ado 
// do disco_plot.ado
// do disco_utils.mata
// do quadprog.ado
// do quadprog.mata


* Step 1: Generate IDs and Time Periods
set obs 20                    // Number of IDs
gen id_col = _n               // Create unique IDs
expand 20                     // Duplicate each ID 20 times (for time periods)
bysort id_col: gen time_col = _n   // Generate time_col within each ID
expand 50                     // Create 50 observations per ID-time pair

* Step 2: Generate group-specific means and standard deviations
bysort id_col time_col: gen double group_mean = runiform()*10 - 5   // Means between -5 and 5
bysort id_col time_col: gen double group_sd   = runiform()*2 + 0.5  // SDs between 0.5 and 2.5

* Step 3: Generate the y_col with group-specific means and variances
gen double y_col = group_mean + group_sd * rnormal()


// Run disco with aggregation
disco y_col id_col time_col, idtarget(1) t0(2) agg("quantile") graph permutation m(100) g(100)
