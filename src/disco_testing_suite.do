clear all
set seed 12345

net install disco, from("C:\Users\david\Downloads\DiSCos_stata\src") replace


//------------------------------------------------------
// Create artificial dataset with varying distributions
//------------------------------------------------------
* Step 1: Generate IDs and Time Periods
set obs 20                    // Number of IDs
gen id = _n                   // Create unique IDs
expand 20                     // Duplicate each ID 20 times (for time periods)
bysort id: gen time = _n      // Generate time within each ID
expand 50                     // Create 50 observations per ID-time pair

* Step 2: Generate group-specific means and standard deviations
bysort id time: gen double group_mean = runiform()*10 - 5   
bysort id time: gen double group_sd = runiform()*2 + 0.5    

* Step 3: Generate the y variable with group-specific means and variances
gen double y = group_mean + group_sd * rnormal()

// Add treatment effect that varies across the distribution
replace y = y + 1 + 0.5*y if id==1 & time>=10  // Treatment effect increases with y


//------------------------------------------------------
// Basic functionality tests
//------------------------------------------------------

// Test 1: Basic disco without options
quietly: disco y id time, idtarget(1) t0(10)

// Check return values
assert e(cmd) == "disco"
assert e(t0) == "10"
assert e(cl) == ".95"
assert e(N) > 0

// Check matrix dimensions
matrix w = e(weights)
mata {
    st_numscalar("dims_ok", cols(st_matrix("w"))==19 & rows(st_matrix("w"))==1)
    assert(st_numscalar("dims_ok"))
}

// Check weights sum to 1 (standard approach with simplex)
mata {
    st_numscalar("sum_w", sum(st_matrix("w")))
    assert(abs(st_numscalar("sum_w") - 1) < 1e-8)
}

// Check non-negativity constraint
mata {
    st_numscalar("min_w", min(st_matrix("w")))
    assert(st_numscalar("min_w") >= -1e-8)
}

//------------------------------------------------------
// Test mixture approach
//------------------------------------------------------

// Test 2: Mixture approach with simplex
quietly: disco y id time, idtarget(1) t0(10) mixture

// Check weights
matrix w_mix = e(weights)
mata {
    st_numscalar("sum_w_mix", sum(st_matrix("w_mix")))
    st_numscalar("min_w_mix", min(st_matrix("w_mix")))
    assert(abs(st_numscalar("sum_w_mix") - 1) < 1e-8)
    assert(st_numscalar("min_w_mix") >= -1e-8)
}

// Check matrix dimensions and content
foreach mat in weights quantile_diff cdf_diff quantile_synth cdf_synth quantile_t cdf_t cids {
    matrix temp = e(`mat')
    mata {
        st_numscalar("has_content", !missing(st_matrix("temp")[1,1]))
        assert(st_numscalar("has_content"))
    }
}

//------------------------------------------------------
// Test without simplex constraint
//------------------------------------------------------

// Test 3: Standard approach without simplex
quietly: disco y id time, idtarget(1) t0(10) nosimplex

// Check weights
matrix w_nosimplex = e(weights)
mata {
    st_numscalar("sum_w_nosimplex", sum(st_matrix("w_nosimplex")))
    assert(abs(st_numscalar("sum_w_nosimplex") - 1) < 1e-8)
}

// Should allow negative weights
mata {
    st_numscalar("min_w_nosimplex", min(st_matrix("w_nosimplex")))
    st_numscalar("has_neg", st_numscalar("min_w_nosimplex") < 0)
    assert(st_numscalar("has_neg") == 1)
}

// Test 4: Mixture approach without simplex
quietly: disco y id time, idtarget(1) t0(10) mixture nosimplex

// Check weights
matrix w_mix_nosimplex = e(weights)
mata {
    st_numscalar("sum_w_mix_nosimplex", sum(st_matrix("w_mix_nosimplex")))
    assert(abs(st_numscalar("sum_w_mix_nosimplex") - 1) < 1e-8)
}

mata {
    st_numscalar("min_w_mix_nosimplex", min(st_matrix("w_mix_nosimplex")))
    st_numscalar("has_neg_mix", st_numscalar("min_w_mix_nosimplex") < 0)
    assert(st_numscalar("has_neg_mix") == 1)
}

//------------------------------------------------------
// Test confidence intervals
//------------------------------------------------------

// Test 5: Confidence intervals
quietly: disco y id time, idtarget(1) t0(10) ci boots(100) cl(0.90)

// Check CI matrices
foreach mat in qdiff_lower qdiff_upper cdiff_lower cdiff_upper {
    matrix temp = e(`mat')
    mata {
        st_numscalar("dims_ok", rows(st_matrix("temp"))==100 & cols(st_matrix("temp"))==20)
        assert(st_numscalar("dims_ok"))
    }
}

//------------------------------------------------------
// Test permutation inference
//------------------------------------------------------

// Test 6: Permutation test
quietly: disco y id time, idtarget(1) t0(10) permutation
assert e(pval) != .
assert e(pval) >= 0 & e(pval) <= 1

//------------------------------------------------------
// Test different grid sizes
//------------------------------------------------------

// Test 7: Different grid sizes
quietly: disco y id time, idtarget(1) t0(10) m(50) g(50)
mata {
    st_numscalar("dims_ok", rows(st_matrix("e(quantile_diff)"))==50 & rows(st_matrix("e(cdf_diff)"))==50)
    assert(st_numscalar("dims_ok"))
}

quietly: disco y id time, idtarget(1) t0(10) m(200) g(200)
mata {
    st_numscalar("dims_ok", rows(st_matrix("e(quantile_diff)"))==200 & rows(st_matrix("e(cdf_diff)"))==200)
    assert(st_numscalar("dims_ok"))
}

//------------------------------------------------------
// Test aggregation and summary statistics
//------------------------------------------------------

// Test 8: Quantile differences
quietly: disco y id time, idtarget(1) t0(10) agg("quantileDiff")
mata {
    st_numscalar("has_content", !missing(st_matrix("e(summary_stats)")))
    assert(st_numscalar("has_content"))
}

// Test 9: CDF differences
quietly: disco y id time, idtarget(1) t0(10) agg("cdfDiff")
mata {
    st_numscalar("has_content", !missing(st_matrix("e(summary_stats)")))
    assert(st_numscalar("has_content"))
}


quietly: disco y id time, idtarget(1) t0(10) agg("cdf") graph

//------------------------------------------------------
// Test disco_plot
//------------------------------------------------------

// Store all matrices for plotting
tempname qd qt qs cd ct cs qdl qdu cdl cdu
matrix `qd' = e(quantile_diff)
matrix `qt' = e(quantile_t)
matrix `qs' = e(quantile_synth)
matrix `cd' = e(cdf_diff)
matrix `ct' = e(cdf_t)
matrix `cs' = e(cdf_synth)
matrix `qdl' = e(qdiff_lower)
matrix `qdu' = e(qdiff_upper)
matrix `cdl' = e(cdiff_lower)
matrix `cdu' = e(cdiff_upper)

// Check matrix dimensions
mata {
    qd = st_matrix(st_local("qd"))
    assert(rows(qd)==100 & cols(qd)==20)
    qt = st_matrix(st_local("qt"))
    assert(rows(qt)==100 & cols(qt)==20)
    qs = st_matrix(st_local("qs"))
    assert(rows(qs)==100 & cols(qs)==20)
}

// Test 11-13: Plot tests
quietly: {
    // Basic quantile plot without CIs
    disco_plot, agg("quantile") m(100) g(100) t_max(20) doci(0) cl(0.95) ///
        quantile_diff(`qd') quantile_t(`qt') quantile_synth(`qs') ///
        cdf_diff(`cd') cdf_synth(`cs') cdf_t(`ct')

    // Quantile differences with CIs
    disco_plot, agg("quantileDiff") m(100) g(100) t_max(20) doci(1) cl(0.95) ///
        quantile_diff(`qd') quantile_t(`qt') quantile_synth(`qs') ///
        cdf_diff(`cd') cdf_synth(`cs') cdf_t(`ct') ///
        qdiff_lower(`qdl') qdiff_upper(`qdu') ///
        cdiff_lower(`cdl') cdiff_upper(`cdu')

    // CDF plot with custom options
    disco_plot, agg("cdf") m(100) g(100) t_max(20) doci(0) cl(0.95) ///
        quantile_diff(`qd') quantile_t(`qt') quantile_synth(`qs') ///
        cdf_diff(`cd') cdf_synth(`cs') cdf_t(`ct') ///
        title("Distribution Effects") ytitle("CDF") xtitle("Y") ///
        color1("blue") color2("red") lwidth("thick")
}

//------------------------------------------------------
// Test grid ranges
//------------------------------------------------------

// Test 14: Check amin and amax are properly stored
quietly: disco y id time, idtarget(1) t0(10)
assert e(amin) < 0  // Given our DGP should be negative
assert e(amax) > 0  // Given our DGP should be positive
assert e(amin) != .
assert e(amax) != .

//------------------------------------------------------
// Test error handling
//------------------------------------------------------

// Tests 15-18: Error handling
rcof "noi disco y id time, idtarget(1) t0(10) agg(invalid)" == 198
rcof "noi disco y id time, idtarget(1) t0(10) ci cl(1.5)" == 198
rcof "noi disco y id time, idtarget(1) t0(10) m(0)" == 198
rcof "noi disco y id time" == 198

//------------------------------------------------------
// Test with missing data
//------------------------------------------------------

// Test 19: Handle missing data
preserve
replace y = . if _n == 1
quietly: disco y id time, idtarget(1) t0(10)
restore

//------------------------------------------------------
// Test different quantile ranges
//------------------------------------------------------

// Test 20: Custom quantile range
quietly: disco y id time, idtarget(1) t0(10) qmin(0.1) qmax(0.9)

display "All tests completed successfully!"
