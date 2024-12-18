capture program drop disco

program define disco, rclass
    version 18.0
    
    // Required: varlist of three variables
    // Required: idtarget and T0
    // Optional: all others
    syntax varlist(min=3 max=3) [if] [in], ///
        idtarget(integer) ///
        T0(integer) ///
        [M(integer 1000) ///
        G(integer 1000) ///
        CI ///
        BOOTS(integer 500) ///
        CL(real 0.95) ///
        QMIN(real 0) ///
        QMAX(real 1) ///
        SIMPLEX ///
        MIXTURE ///
        PERMUTATION ///
        SEED(integer -1) ///
        UNIFORM]
        
    // Mark the estimation sample
    marksample touse, novarlist
    
    // Extract variables from varlist
    local y_col : word 1 of `varlist'
    local id_col : word 2 of `varlist'
    local time_col : word 3 of `varlist'
    
    // Mark out missing values in required variables
    markout `touse' `y_col' `time_col'
    markout `touse' `id_col', strok  // strok option for string IDs if needed
    
    // Initialize optional arguments with defaults if not specified
    if ("`m'"=="") local m = 1000
    if ("`g'"=="") local g = 1000
    if ("`boots'"=="") local boots = 500
    if ("`cl'"=="") local cl = 0.95
    if ("`qmin'"=="") local qmin = 0
    if ("`qmax'"=="") local qmax = 1
    if ("`seed'"=="") local seed = -1
    
    // Check required numeric options
    if missing(`t0') {
        di as err "t0() is required"
        exit 198
    }
    if missing(`idtarget') {
        di as err "idtarget() is required"
        exit 198
    }
    

    // Initialize all flags to 0 by default
    local simplex_flag = 0
    local mixture_flag = 0
    local permutation_flag = 0
    local doci = 0
    local uniform_flag = 0

    // Only set flags to 1 if the corresponding option was specified
    if "`simplex'" != "" {
        local simplex_flag = 1
    }
    if "`mixture'" != "" {
        local mixture_flag = 1
    }
    if "`permutation'" != "" {
        local permutation_flag = 1
    }
    if "`ci'" != "" {
        local doci = 1
    }
    if "`uniform'" != "" {
        local uniform_flag = 1
    }

    // Extract time range
    quietly levelsof `time_col', local(times)
    local min_time : word 1 of `times'
    local max_time : word `=wordcount("`times'")' of `times'
    gen t_col = `time_col' - `min_time' + 1
    local t_max = `max_time' - `min_time' + 1
    local t0_col = `t0' - `min_time' + 1
	


    // Argument checks
    if `m' < 1 {
        di as err "M must be >=1"
        exit 198
    }
    if `g' < 2 {
        di as err "G must be >=2"
        exit 198
    }
    if `qmin' < 0 | `qmax' > 1 {
        di as err "q_min must be >=0 and q_max <=1"
        exit 198
    }

    // Preserve dataset before modifications
    tempname base
    preserve
    
    // Keep only the estimation sample
    keep if `touse'

    mata: mata clear
    // Include the MATA utility file that contains all the functions
    // Make sure disco_utils.mata is accessible in your adopath or same folder
    include disco_utils.mata
	include quadprog.ado
	include quadprog.mata

    mata {
        M = `m'
        G = `g'
        T0 = `t0_col'
        T_max = `t_max'
        q_min = `qmin'
        q_max = `qmax'
        cl = `cl'
        nboots = `boots'
        simplex = `simplex_flag'
        mixture = `mixture_flag'
        uniform = `uniform_flag'

        y = st_data(.,"`y_col'")
        id = st_data(.,"`id_col'")
        tt = st_data(.,"t_col")
        target_id = `idtarget'

        // Run the main procedure 
        // Execute the main analysis in Mata
        rc = disco_wrapper( ///
            y, id, tt, ///
            target_id, T0, T_max, M, G, ///
            q_min, q_max, simplex, mixture)
            
        if (`permutation_flag'==1) {
            pval = disco_permutation_test(y,id,tt,target_id,T0,T_max,M,G,q_min,q_max,simplex,mixture)
            st_local("pval", strofreal(pval))
        }
        
        if (`doci'==1) {
            rc2 = disco_ci_wrapper(y, id, tt, target_id, T0, T_max, M, G,
                                q_min, q_max, simplex, mixture,
                                nboots, cl, uniform)
            st_local("rc2", strofreal(rc2))

// 			if (`rc2' != 0) {
// 				error `rc2'
// 			}
        }

    }
    
//     if `rc' != 0 {
//         error `rc'
//     }
    restore
	
	di("`doci'")
	
	// return function output
	if `doci' == 1 {
		return matrix qdiff_lower = qdiff_lower
		return matrix qdiff_upper = qdiff_upper
		return matrix cdiff_lower = cdiff_lower
		return matrix cdiff_upper = cdiff_upper
	}
	return matrix weights = weights
	return matrix quantile_diff = quantile_diff
	return matrix cdf_diff = cdf_diff

	
    
end
