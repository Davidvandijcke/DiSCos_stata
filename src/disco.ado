capture program drop disco

program define disco, rclass
    version 18.0
    

    // Required: varlist of three variables
    // Required: idtarget and T0
    // Optional: all others
    syntax varlist(min=3 max=3) [if] [in], ///
        idtarget(integer) ///
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

    // Extract variables from varlist
    local y_col : word 1 of `varlist'
    local id_col : word 2 of `varlist'
    local time_col : word 3 of `varlist'

    // Check required numeric options
    if missing(`T0') {
        di as err "T0() is required"
        exit 198
    }
    if missing(`idtarget') {
        di as err "idtarget() is required"
        exit 198
    }

    // Set defaults for optional arguments if not specified
    if missing(`M') local M = 1000
    if missing(`G') local G = 1000
    if missing(`BOOTS') local BOOTS = 500
    if missing(`CL') local CL = 0.95
    if missing(`QMIN') local QMIN = 0
    if missing(`QMAX') local QMAX = 1

    // For CI, SIMPLEX, MIXTURE, PERMUTATION, UNIFORM:
    // If the option is not specified, the corresponding local is empty.
    // Assign default values based on presence/absence of the option.
    if "`CI'" == "" local CI "no"
    else local CI "yes"

    if "`SIMPLEX'" == "" local SIMPLEX "no"
    else local SIMPLEX "yes"

    if "`MIXTURE'" == "" local MIXTURE "no"
    else local MIXTURE "yes"

    // PERMUTATION is a toggle, if specified, it's on; if not, off.
    local permutation_flag = 0
    if "`PERMUTATION'" != "" local permutation_flag = 1

    if "`UNIFORM'" == "" local UNIFORM "no"
    else local UNIFORM "yes"

    // Set the seed if provided
    if `SEED' != . {
        set seed `SEED'
    }
		 
// 	// TODO: for testing, comment out
// 	clear
// 	set seed 123
// 	set obs 20       // Number of IDs
// 	gen id_col = _n
// 	expand 20        // Duplicate each ID 20 times (for time periods)
// 	bysort id: gen time_col = _n
// 	expand 50        // Create 50 observations per ID-time pair
// 	gen y_col = rnormal()
//
//		
//	
// 	local varlist y_col id_col time_col
// 	local idtarget = 1
// 	local T0 = 2
// 	local M = 20
// 	local G = 10
// 	local CI = "yes"
// 	local boots = 10
// 	local cl = 0.95
// 	local qmin = 0
// 	local qmax=1
// 	local SEED = 124
// 	local mixture = "no"
// 	local permutation = "no"
// 	local uniform = "no"

    // Extract variable names from varlist
    local y_col : word 1 of `varlist'
    local id_col : word 2 of `varlist'
    local time_col : word 3 of `varlist'

    // Check mandatory arguments
    if missing(`T0') {
        di as err "T0() is required"
        exit 198
    }
    if missing(`idtarget') {
        di as err "idtarget() is required"
        exit 198
    }

    // Check if seed is provided
    if `SEED'!=. {
        set seed `SEED'
    }

    // Extract time range
    quietly levelsof `time_col', local(times)
    local min_time : word 1 of `times'
    local max_time : word `=wordcount("`times'")' of `times'
    gen t_col = `time_col' - `min_time' + 1
    local T_max = `max_time' - `min_time' + 1
    local T0_col = `T0' - `min_time' + 1

    // Argument checks
    if `M' < 1 {
        di as err "M must be >=1"
        exit 198
    }
    if `G' < 2 {
        di as err "G must be >=2"
        exit 198
    }
    if `qmin' < 0 | `qmax' > 1 {
        di as err "q_min must be >=0 and q_max <=1"
        exit 198
    }

    // Convert yes/no strings to numeric flags

    local simplex_flag = 0
    if "`simplex'"=="yes" {
        local simplex_flag = 1
    }

    local mixture_flag = 0
    if "`mixture'"=="yes" {
        local mixture_flag = 1
    }

    local perm_flag = 0
    if "`permutation'"=="yes" {
        local perm_flag = 1
    }
	
	local doCI = 0
	if "`CI'"=="yes" {
		local doCI = 1
	}

    local uniform_flag = 0
    if "`uniform'"=="yes" {
        local uniform_flag = 1
    }

    // Preserve dataset
    tempname base
    preserve
    save `base', replace

    mata: mata clear
    // Include the MATA utility file that contains all the functions
    // Make sure disco_utils.mata is accessible in your adopath or same folder
    include disco_utils.mata

    mata: 
        M = `M'
        G = `G'
        T0 = `T0_col'
        T_max = `T_max'
        q_min = `qmin'
        q_max = `qmax'
        cl = `CL'
        nboots = `boots'
        simplex = `simplex_flag'
        mixture = `mixture_flag'
        uniform = `uniform_flag'

        y = st_data(.,"`y_col'")
        id = st_data(.,"`id_col'")
        tt = st_data(.,"t_col")
        target_id = `idtarget'

        // Run the main procedure
        results = disco_full_run(y,id,tt,target_id,T0,T_max,M,G,q_min,q_max,simplex,mixture)

        // Store main outputs
        st_matrix("r(weights)", results.weights')
        st_matrix("r(quantile_diff)", results.quantile_diff)
        st_matrix("r(cdf_diff)", results.cdf_diff)
	end
	
	if "`permutation'"=="yes" {
		mata: 
			pval = disco_permutation_test(y,id,tt,target_id,T0,T_max,M,G,q_min,q_max,simplex,mixture)
			st_local("pval", strofreal(pval))
		end
	}
	
	if "`CI'"=="yes" {
		mata: 
			CIres = disco_bootstrap_CI(y,id,tt,target_id,T0,T_max,M,G,q_min,q_max,simplex,mixture,nboots,cl,uniform)
			st_matrix("r(qdiff_lower)", CIres.qdiff_lower)
			st_matrix("r(qdiff_upper)", CIres.qdiff_upper)
			st_matrix("r(cdiff_lower)", CIres.cdiff_lower)
			st_matrix("r(cdiff_upper)", CIres.cdiff_upper)	
		end	
	}

}
