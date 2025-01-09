/*
Distributional Synthetic Controls (DiSCo)
Implementation based on Gunsilius (2023) and Van Dijcke, Gunsilius, and Wright (2024)

This command implements the DiSCo method for estimating counterfactual distributions
using synthetic controls. It allows for both quantile-based and mixture-based approaches,
with options for confidence intervals, permutation tests, and summary statistics.

Syntax:
    disco varlist(3) [if] [in], idtarget(integer) t0(integer) [options]

Required input:
    varlist:    Three variables in order:
                1. Outcome variable (numeric)
                2. Unit ID variable (numeric)
                3. Time period variable (integer)
    idtarget:   ID of the treated unit
    t0:         First treatment period

Options:
	M(integer): 		Number of quantile samples to use in estimation for approximating the integral (default: 1000).
    G(integer):         Number of grid points to evaluate quantile function/cdf (default: 100)
    CI:                 Compute confidence intervals
    BOOTS(integer):     Number of bootstrap replications (default: 300)
    CL(real):          Confidence level (default: 0.95)
    QMIN(real):        Minimum quantile for estimation (default: 0)
    QMAX(real):        Maximum quantile for estimation (default: 1)
    NOSIMPLEX:         Don't constrain weights to unit simplex
    MIXTURE:           Use mixture of distributions approach
    PERMUTATION:       Perform permutation test
    SEED(integer):     Set random seed
    NOUNIFORM:         Don't use uniform confidence bands
    AGG(string):       Type of aggregation for summary statistics
                       ("quantile", "cdf", "quantileDiff", "cdfDiff"). 
    SAMPLES(numlist):  Quantile points for summary statistics

Stored results:
    e(weights):        Synthetic control weights
    e(quantile_diff):  Quantile differences
    e(cdf_diff):       CDF differences
    e(quantile_synth): Synthetic control quantiles
    e(quantile_t):     Target unit quantiles
    e(cdf_synth):      Synthetic control CDFs
    e(cdf_t):          Target unit CDFs
    e(summary_stats):  Summary statistics (if agg specified)
    
    If CI specified:
    e(qdiff_lower):    Lower bound for quantile differences
    e(qdiff_upper):    Upper bound for quantile differences
    e(cdiff_lower):    Lower bound for CDF differences
    e(cdiff_upper):    Upper bound for CDF differences

Author: David Van Dijcke
Version: 1.0.0
Date: December 2024
*/


program define disco, eclass
    version 18.0
	

    
    // Syntax parsing
    syntax varlist(min=3 max=3) [if] [in], ///
        idtarget(integer) ///
        T0(integer) ///
        [M(integer 100) ///
        G(integer 100) ///
        CI ///
        BOOTS(integer 300) ///
        CL(real 0.95) ///
        QMIN(real 0) ///
        QMAX(real 1) ///
        NOSIMPlex ///
        MIXture ///
        PERMutation ///
        SEED(integer -1) ///
        NOUNIForm ///
        AGG(string) ///
        SAMPles(numlist)]
    
	
    // Input validation
    if !inlist("`agg'", "", "quantile", "cdf", "quantileDiff", "cdfDiff") {
        di as error "agg() must be one of: quantile, cdf, quantileDiff, cdfDiff"
        exit 198
    }
	
// 	// load mata objects
// 	mata mata mlib index
    
    // Mark the estimation sample
    marksample touse, novarlist
    
    // Extract variable names from varlist
    local y_col : word 1 of `varlist'
    local id_col : word 2 of `varlist'
    local time_col : word 3 of `varlist'
    
    // Check for missing values
    markout `touse' `y_col' `time_col'
    markout `touse' `id_col', strok
    
    // Initialize optional arguments
	if ("`m'"=="") local m = 100
    if ("`g'"=="") local g = 100
    if ("`boots'"=="") local boots = 300
    if ("`cl'"=="") local cl = 0.95
    if ("`qmin'"=="") local qmin = 0
    if ("`qmax'"=="") local qmax = 1
    if ("`seed'"=="") local seed = -1
    if ("`samples'"=="") local samples "0.25 0.5 0.75"
    
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
    local simplex_flag = 1
    local mixture_flag = 0
    local permutation_flag = 0
    local doci = 0
    local uniform_flag = 1

    if "`nosimplex'" != "" local simplex_flag = 0
    if "`mixture'" != "" local mixture_flag = 1
    if "`permutation'" != "" local permutation_flag = 1
    if "`ci'" != "" local doci = 1
    if "`nouniform'" != "" local uniform_flag = 0
	
    
    // Additional validation checks
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
	if `cl' < 0 | `cl' > 1 {
		di as err "cl must be >=0 and <=1"
        exit 198
	}
	// Preserve dataset before Mata operations
    tempname base
    preserve
    quietly: keep if `touse'
    
    // Identify time range in data
    quietly levelsof `time_col', local(times)
    local min_time : word 1 of `times'
    local max_time : word `=wordcount("`times'")' of `times'
    gen t_col = `time_col' - `min_time' + 1
    local t_max = `max_time' - `min_time' + 1
    local t0_col = `t0' - `min_time' + 1
    

    
    //************************
    // Main analysis in Mata
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
        
        // Load data into Mata
        y = st_data(.,"`y_col'")
        id = st_data(.,"`id_col'")
        tt = st_data(.,"t_col")
        target_id = `idtarget'
        
		
        // Run main DiSCo analysis
			
        rc = disco_wrapper(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)
        
        // Permutation test if requested
        if (`permutation_flag'==1) {
            pval = disco_permutation_test(y,id,tt,target_id,T0,T_max, M, G,q_min,q_max,simplex,mixture)
            st_local("pval", strofreal(pval))
        };
        
        // Confidence intervals if requested
        if (`doci'==1) {
            rc2 = disco_ci_wrapper(y, id, tt, target_id, T0, T_max, M, G,
                                q_min, q_max, simplex, mixture,
                                nboots, cl, uniform, st_matrix("quantile_diff"), st_matrix("cdf_diff"))
            st_local("rc2", strofreal(rc2))
        };
        
        // Summary statistics if requested
		if (anyof(("quantile", "cdf"), "`agg'")) {
				printf("Levels requested so no summary stats table produced.")
		}
        else if ("`agg'" != "") {

            samples_str = "`samples'"
            sample_points = strtoreal(tokens(samples_str))
			quantile_diff_mata = st_matrix("quantile_diff")
			cdf_diff_mata = st_matrix("cdf_diff")
			           
            rc3 = compute_summary_stats("`agg'", sample_points, T0, T_max, quantile_diff_mata,
			cdf_diff_mata, `doci', cl)
        };
    }
	//************************

//	
// 	// Generate plots if requested
//     if "`graph'" != "" & "`agg'" != "" {
//         tempname qd qt qs cd cs qdl qdu cdl cdu 
//         matrix `qd' = quantile_diff
//         matrix `qt' = quantile_t
//         matrix `qs' = quantile_synth
//         matrix `cd' = cdf_diff
//         matrix `cs' = cdf_synth
//
// 		local amin = amin
// 		local amax = amax
//        
//         if `doci' == 1 {
//             matrix `qdl' = qdiff_lower
//             matrix `qdu' = qdiff_upper
//             matrix `cdl' = cdiff_lower
//             matrix `cdu' = cdiff_upper
//         }
//         quietly: disco_plot, agg("`agg'") m(`m') g(`g') t_max(`t_max') doci(`doci') cl(`cl') ///
//             quantile_diff(`qd') quantile_t(`qt') quantile_synth(`qs') ///
//             cdf_diff(`cd') cdf_synth(`cs') cdf_t(cdf_t) ///
//             qdiff_lower(`qdl') qdiff_upper(`qdu') cdiff_lower(`cdl') cdiff_upper(`cdu') ///
//             xmin(`amin') xmax(`amax') `options'
//     }
//	
    // Store results
    if `doci' == 1 {
        ereturn matrix qdiff_lower = qdiff_lower
        ereturn matrix qdiff_upper = qdiff_upper
        ereturn matrix cdiff_lower = cdiff_lower
        ereturn matrix cdiff_upper = cdiff_upper
    }
    
    ereturn matrix quantile_diff = quantile_diff
    ereturn matrix cdf_diff = cdf_diff
    ereturn matrix quantile_synth = quantile_synth
    ereturn matrix quantile_t = quantile_t
    ereturn matrix cdf_synth = cdf_synth
    ereturn matrix cdf_t = cdf_t
	ereturn scalar amin = amin
	ereturn scalar amax = amax

	if !inlist("`agg'", "quantile", "cdf") & "`agg'" != "" {
		ereturn matrix summary_stats = summary_stats
	}

    
    // Store metadata
    ereturn local cmd "disco"
    ereturn local cmdline `"disco `0'"'
    ereturn local agg "`agg'"
    ereturn local cl = `cl'
    ereturn local t0 = `t0'
	ereturn scalar m = `m'
	ereturn scalar g = `g'
	ereturn scalar t_max = `t_max'
	ereturn local doci = `doci'
    ereturn scalar N = _N
    

    // Display permutation test results if requested
    if `permutation_flag' {
        di _n as txt "Permutation test p-value: " as res %5.3f `pval'
		ereturn scalar pval = `pval'
    } 
	else {
		ereturn scalar pval = .
	}
	
	ereturn matrix weights = weights
	ereturn matrix cids = cids

	

end
