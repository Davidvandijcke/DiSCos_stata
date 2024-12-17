// Test data section
clear
set seed 123
set obs 20       // Number of IDs
gen id_col = _n
expand 20        // Duplicate each ID 20 times (for time periods)
bysort id: gen time_col = _n
expand 50        // Create 50 observations per ID-time pair
gen y_col = rnormal()

// Test locals
local varlist y_col id_col time_col
local idtarget = 1
local T0 = 2
local M = 20
local G = 10
local CI = "yes"
local boots = 10
local cl = 0.95
local qmin = 0
local qmax=1
local SEED = 124
local mixture = "no"
local permutation = "no"
local uniform = "no"

capture program drop disco

program define disco, rclass
    version 18.0
    
    syntax varlist(min=3 max=3) [if] [in] ///
        , idtarget(integer) ///
        T0(integer) ///
        [M(integer)] ///
        [G(integer)] ///
        [CI] ///
        [BOOTS(integer)] ///
        [CL(real)] ///
        [QMIN(real)] ///
        [QMAX(real)] ///
        [SIMPLEX] ///
        [MIXTURE] ///
        [PERMUTATION] ///
        [SEED(integer)] ///
        [UNIFORM]

	di("hello")

end

// Now call the program with test data
disco y_col id_col time_col
