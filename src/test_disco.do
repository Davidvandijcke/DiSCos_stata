clear
set seed 123
set obs 20       // Number of IDs
gen id_col = _n
expand 20        // Duplicate each ID 20 times (for time periods)
bysort id_col: gen time_col = _n
expand 50        // Create 50 observations per ID-time pair
gen y_col = rnormal()



// Test full CI with corrected syntax
disco y_col id_col time_col, idtarget(1) T0(2) M(20) G(10) CI boots(10) ///
    cl(0.90) qmin(0) qmax(1) simplex mixture permutation ///
    uniform seed(123)

matrix list r(qdiff_lower)

