clear
set seed 123
set obs 20       // Number of IDs
gen id_col = _n
expand 20        // Duplicate each ID 20 times (for time periods)
bysort id_col: gen time_col = _n
expand 50        // Create 50 observations per ID-time pair
gen y_col = rnormal()

capture program drop disco
include disco.ado


// Test full CI with corrected syntax
disco y_col id_col time_col, idtarget(1) t0(2) simplex
return list
matrix list r(weights)
