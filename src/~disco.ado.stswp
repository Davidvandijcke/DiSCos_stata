capture program drop disco

program define disco, rclass
    version 18.0
    
    // Syntax parsing
    // varlist must have 3 variables: outcome, id, time
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
    
    // Extract variable names from varlist
    local y_col : word 1 of `varlist'
    local id_col : word 2 of `varlist'
    local time_col : word 3 of `varlist'
    
    // Check for missing values
    markout `touse' `y_col' `time_col'
    markout `touse' `id_col', strok
    
    // Initialize optional arguments if not specified
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

    // Set flags based on options
    local simplex_flag = 0
    local mixture_flag = 0
    local permutation_flag = 0
    local doci = 0
    local uniform_flag = 0

    if "`simplex'" != "" local simplex_flag = 1
    if "`mixture'" != "" local mixture_flag = 1
    if "`permutation'" != "" local permutation_flag = 1
    if "`ci'" != "" local doci = 1
    if "`uniform'" != "" local uniform_flag = 1

    // Identify time range in data
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

    // Preserve dataset before Mata operations
    tempname base
    preserve
    
    keep if `touse'

    mata: mata clear
    // Include the Mata file with utility functions and QP solver
    include disco_utils.mata
    include quadprog.ado
    include quadprog.mata

    mata {
        // Store options in Mata variables
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

        // Run main DiSCo analysis
        rc = disco_wrapper(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)
            
        // If permutation test requested
        if (`permutation_flag'==1) {
            pval = disco_permutation_test(y,id,tt,target_id,T0,T_max,M,G,q_min,q_max,simplex,mixture)
            st_local("pval", strofreal(pval))
        }
        
        // If confidence intervals requested
        if (`doci'==1) {
            rc2 = disco_ci_wrapper(y, id, tt, target_id, T0, T_max, M, G,
                                q_min, q_max, simplex, mixture,
                                nboots, cl, uniform)
            st_local("rc2", strofreal(rc2))
        }

    }
    restore

    // If CI was computed, return intervals
    if `doci' == 1 {
        return matrix qdiff_lower = qdiff_lower
        return matrix qdiff_upper = qdiff_upper
        return matrix cdiff_lower = cdiff_lower
        return matrix cdiff_upper = cdiff_upper
        
    }

    // If no CI requested, just return main results
    return matrix weights = weights
    return matrix quantile_diff = quantile_diff
    return matrix cdf_diff = cdf_diff
    return matrix quantile_synth = quantile_synth
    return matrix quantile_t = quantile_t
    return matrix cdf_synth = cdf_synth
    return matrix cdf_t = cdf_t
end
