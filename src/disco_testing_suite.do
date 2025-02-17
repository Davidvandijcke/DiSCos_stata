************************************************************
* Install disco (only once, outside the tests)
************************************************************
clear all
mata: mata clear
cd "/Users/davidvandijcke/University of Michigan Dropbox/David Van Dijcke/Flo_GSRA/stata_repo/src"

net install disco, from("/Users/davidvandijcke/University of Michigan Dropbox/David Van Dijcke/Flo_GSRA/stata_repo/src") replace

// net install disco, from("https://raw.githubusercontent.com/Davidvandijcke/DiSCos_stata/dev/src/") replace

************************************************************
* Now we run each test, calling gen_data prior to each test
************************************************************

*----------------------------------------------------------------------
* Test 1: Basic disco without options
*----------------------------------------------------------------------
gen_data
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


*----------------------------------------------------------------------
* Test 2: Mixture approach with simplex
*----------------------------------------------------------------------
gen_data
disco y id time, idtarget(1) t0(10) mixture

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


*----------------------------------------------------------------------
* Test 3: Standard approach without simplex
*----------------------------------------------------------------------
gen_data
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


*----------------------------------------------------------------------
* Test 4: Mixture approach without simplex
*----------------------------------------------------------------------
gen_data
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


*----------------------------------------------------------------------
* Test 5: Confidence intervals
*----------------------------------------------------------------------
gen_data
// mata: mata clear
// do disco_utils.mata
disco y id time, idtarget(1) t0(10) ci boots(34) cl(0.90) 
disco_plot

// Check CI matrices
foreach mat in qdiff_lower qdiff_upper cdiff_lower cdiff_upper {
    matrix temp = e(`mat')
    mata {
        st_numscalar("dims_ok", rows(st_matrix("temp"))==100 & cols(st_matrix("temp"))==20)
        assert(st_numscalar("dims_ok"))
    }
}


*----------------------------------------------------------------------
* Test 6: Permutation test
*----------------------------------------------------------------------
gen_data
quietly: disco y id time, idtarget(1) t0(10) permutation
assert e(pval) != .
assert e(pval) >= 0 & e(pval) <= 1


*----------------------------------------------------------------------
* Test 7: Different grid sizes
*----------------------------------------------------------------------
gen_data
quietly: disco y id time, idtarget(1) t0(10)  g(50)
mata {
    st_numscalar("dims_ok", rows(st_matrix("e(quantile_diff)"))==50 & rows(st_matrix("e(cdf_diff)"))==50)
    assert(st_numscalar("dims_ok"))
}

gen_data
quietly: disco y id time, idtarget(1) t0(10) g(100) ci boots(100)
mata {
    st_numscalar("dims_ok", rows(st_matrix("e(quantile_diff)"))==100 & rows(st_matrix("e(cdf_diff)"))==100)
    assert(st_numscalar("dims_ok"))
}


*----------------------------------------------------------------------
* Test 8: Aggregation "quantileDiff"
*----------------------------------------------------------------------

gen_data
quietly: disco y id time, idtarget(1) t0(10) agg("quantileDiff")  
mata {
    st_numscalar("has_content", !(eltype(st_matrix("e(summary_stats)")) == ""))
    assert(st_numscalar("has_content"))
}


*----------------------------------------------------------------------
* Test 9: CDF differences
*----------------------------------------------------------------------
gen_data
disco y id time, idtarget(1) t0(10) agg("cdfDiff") 
mata {
    st_numscalar("has_content", !(eltype(st_matrix("e(summary_stats)")) == "" ))
    assert(st_numscalar("has_content"))
}


gen_data
disco y id time, idtarget(1) t0(10) agg("cdf") 
// TODO: fix - local gmin = `min' = local gmin =  ...

*----------------------------------------------------------------------
* Tests 10-13: Disco_plot
*----------------------------------------------------------------------

gen_data
quietly: disco y id time, idtarget(1) t0(10) ci boots(100) // so we have some CI mats stored


// Test 10: CDF plot with custom options
quietly: disco_plot, title("Distribution Effects") ytitle("CDF") xtitle("Y") ///
    color1("blue") color2("red") lwidth("thick")
	
	
// with cdfDiff option
gen_data

// Test 11: CDF diff plot with custom options
quietly: disco y id time, idtarget(1) t0(10) agg("cdfDiff") // so we have some CI mats stored
disco_plot

// Test 12: CDF plot
quietly: disco y id time, idtarget(1) t0(10) agg("cdf") // so we have some CI mats stored
disco_plot

// Test 13: quantile plot
quietly: disco y id time, idtarget(1) t0(10) agg("quantile") // so we have some CI mats stored
disco_plot

// Test 13.1 different plot than specified
gen_data
quietly: disco y id time, idtarget(1) t0(10) agg("quantile") // so we have some CI mats stored
disco_plot, agg("quantileDiff")



*----------------------------------------------------------------------
* Test 14: Check amin and amax
*----------------------------------------------------------------------
gen_data
quietly: disco y id time, idtarget(1) t0(10)
assert e(amin) < 0
assert e(amax) > 0
assert e(amin) != .
assert e(amax) != .


*----------------------------------------------------------------------
* Tests 15-18: Error handling
*----------------------------------------------------------------------
gen_data
rcof "noi disco y id time, idtarget(1) t0(10) agg(invalid)" == 198
rcof "noi disco y id time, idtarget(1) t0(10) ci cl(1.5)" == 198
rcof "noi disco y id time, idtarget(1) t0(10) g(0)" == 198
rcof "noi disco y id time" == 198


*----------------------------------------------------------------------
* Test 19: Missing data
*----------------------------------------------------------------------
gen_data
preserve
replace y = . if _n == 1
quietly: disco y id time, idtarget(1) t0(10)
restore


*----------------------------------------------------------------------
* Test 20: Custom quantile range
*----------------------------------------------------------------------
gen_data
quietly: disco y id time, idtarget(1) t0(10) qmin(0.1) qmax(0.9)

display "All tests completed successfully!"

*----------------------------------------------------------------------
* Test 21: Post-estimation table
*----------------------------------------------------------------------
gen_data
quietly: disco y id time, idtarget(1) t0(10) agg("quantileDiff") 
disco_estat summary


*----------------------------------------------------------------------
* Test 22: weight function
*----------------------------------------------------------------------
gen_data
gen str_id = "control " + string(id)
disco y id time, idtarget(1) t0(3) ci boots(100) cl(0.95) agg("quantileDiff")
disco_estat summary
disco_plot
disco_weight id str_id
disco_weight id str_id, n(10) round(0.01)




