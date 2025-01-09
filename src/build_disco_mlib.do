clear all
cd "/Users/davidvandijcke/University of Michigan Dropbox/David Van Dijcke/Flo_GSRA/stata_repo/src"

mata:
mata clear
mata set matastrict off
end

// Load all Mata functions into memory
do disco_utils.mata
do quadprog.mata

// Now create and populate the Mata library
mata:
mata mlib create ldisco, replace

// Add each function you defined in disco_utils.mata
mata mlib add ldisco CI_out()
mata mlib add ldisco iter_out()
mata mlib add ldisco boot_out()
mata mlib add ldisco disco_out()
mata mlib add ldisco get_unique()
mata mlib add ldisco disco_quantile_points()
mata mlib add ldisco disco_quantile()
mata mlib add ldisco disco_solve_weights()
mata mlib add ldisco disco_mixture_weights()
mata mlib add ldisco cdf_builder()
mata mlib add ldisco cdf_at_points()
mata mlib add ldisco disco_full_run()
mata mlib add ldisco disco_compute_ratio()
mata mlib add ldisco disco_permutation_test()
mata mlib add ldisco disco_CI_iter()
mata mlib add ldisco bootCounterfactuals()
mata mlib add ldisco disco_bootstrap_CI()
mata mlib add ldisco disco_wrapper()
mata mlib add ldisco disco_ci_wrapper()
mata mlib add ldisco compute_summary_stats()
mata mlib add ldisco solve_quadprog()
end
//
// // Done. The library ldisco.mlib is now created with all functions.
//
//
// // Output structure for confidence intervals
// struct CI_out {
//     real matrix qdiff_lower,    // Lower bound for quantile differences
//               qdiff_upper,      // Upper bound for quantile differences
//               cdiff_lower,      // Lower bound for CDF differences
//               cdiff_upper       // Upper bound for CDF differences
// }
//
// // Output structure for bootstrap iteration
// struct iter_out {
//     real vector target_q,       // Target quantiles
//               target_cdf,       // Target CDF
//               weights          // Optimal weights
//     real matrix controls_q,     // Control unit quantiles
//               controls_cdf      // Control unit CDFs
// }
//
// // Output structure for bootstrap results
// struct boot_out {
//     real matrix quantile_diff,  // Quantile differences
//               cdf_diff,         // CDF differences
//               quantile_synth,   // Synthetic quantiles
//               cdf_synth,        // Synthetic CDFs
//               quantile_t,       // Target quantiles
//               cdf_t             // Target CDFs
// }
//
// struct disco_out {
//     real matrix weights, quantile_diff, cdf_diff, quantile_synth, cdf_synth,  
//     quantile_t, cdf_t, cids
// }
