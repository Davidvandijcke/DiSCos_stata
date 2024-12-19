clear
set seed 123

capture program drop _all
mata: mata clear

include disco_utils.mata
include quadprog.ado
include quadprog.mata
include disco_estat.ado
include disco.ado
include disco_plot.ado


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
