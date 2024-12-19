mata:


// Returns unique values from a vector
real vector unique(real vector x) {
    real vector y
    y = uniqrows(sort(x,1))
    return(y)
}

// Compute quantiles at arbitrary probabilities p using type=7 interpolation
// Parameters:
//   X: Input vector of values
//   p: Vector of probability points at which to compute quantiles
real vector disco_quantile_points(real vector X, real vector p) {
    real vector Xs, out
    real scalar N, prob, alpha, floor_alpha, gamma, idx
    
    N = length(X)
    if (N==0) {
        out = J(length(p),1,.)
        return(out)
    }

    Xs = sort(X,1)
    out = J(length(p),1,.)
    
    // Calculate quantile for each probability point
    for (i=1; i<=length(p); i++) {
        prob = p[i]
        if (prob<=0) {
            out[i] = Xs[1]
        }
        else if (prob>=1) {
            out[i] = Xs[N]
        } else {
            alpha = (N-1)*prob + 1
            floor_alpha = floor(alpha)
            gamma = alpha - floor_alpha
            
            if (floor_alpha<1) {
                out[i] = Xs[1]
            } else if (floor_alpha>=N) {
                out[i] = Xs[N]
            } else {
                idx = floor_alpha
                out[i] = Xs[idx]*(1-gamma) + Xs[idx+1]*gamma
            }
        }
    }
    return(out)
}

// Compute quantiles at M evenly spaced points between q_min and q_max
// Parameters:
//   X: Input vector
//   M: Number of quantile points
//   q_min, q_max: Range for quantile computation
real vector disco_quantile(real vector X, real scalar M, real scalar q_min, real scalar q_max) {
    real vector p
    real scalar j, p_j
    
    // Generate evenly spaced probability points
    p = J(M,1,.)
    for (j=1; j<=M; j++) {
        p_j = q_min + (q_max - q_min)*(j-1)/(M-1)
        p[j] = p_j
    }
    
    return(disco_quantile_points(X,p))
}

// Solve for optimal weights using quadratic programming
// Parameters:
//   controls: Matrix of control variables
//   target: Target vector to match
//   M: Number of quantile points
//   q_min, q_max: Range for quantile computation
//   simplex: Whether to constrain weights to be non-negative
real vector disco_solve_weights(real matrix controls,
                              real vector target,
                              real scalar M, q_min, q_max,
                              string scalar qmethod,
                              real scalar qtype, simplex)
{
    real scalar J, i, sc
    real matrix controls_s, C_in, G, CE, CI
    real vector target_s, p, g0, ce0, ci0, w, res
    
    J = cols(controls)
    
    // Construct probability vector
    p = J(M,1,.)
    for (i=1; i<=M; i++) {
        p[i] = q_min + (q_max - q_min)*(i-1)/(M-1)
    }

    // Compute quantiles for controls and target
    real matrix Cq
    Cq = J(M,J,.)
    for (i=1; i<=J; i++) {
        Cq[,i] = disco_quantile_points(controls[,i], p)
    }
    target_s = disco_quantile_points(target, p)

    // Setup QP problem matrices
    C_in = Cq  // No scaling applied
    G = 2*(C_in' * C_in)    // J × J matrix
    g0 = -2*(C_in' * target_s)  // J × 1 vector
    
    // Equality constraint: sum of weights = 1
    CE = J(J,1,1)     // J × 1 matrix
    ce0 = (-1)        // 1 × 1 vector
    
    // Inequality constraints
    if (simplex==1) {
        CI = I(J)         // J × J matrix for non-negativity constraints
        ci0 = J(J,1,0)    // J × 1 vector of zeros
    } else {
        // Single non-binding constraint: sum of weights >= -inf -- because of how I've set up the quadprog plugin
        // This effectively imposes no constraint since we already have sum of weights = 1
        CI = J(J,1,1)     // J × 1 matrix of ones (for sum of weights)
        ci0 = (-1e20)     // 1 × 1 vector with large negative number
    }

    // Solve QP and extract weights
    res = solve_quadprog(G, g0, CE, ce0, CI, ci0)
    w = res[1..J]

    return(w)
}
// Compute mixture weights using linear programming
// Parameters:
//   controls: Matrix of control variables
//   target: Target vector to match
//   M: Number of mixture points
//   simplex: Whether to constrain weights to be non-negative
real vector disco_mixture_weights(real matrix controls, real vector target, real scalar M, simplex) {
    // Declare all scalars
    real scalar i, J, n, grid_min, grid_max, ss, val, g, total_vars
    
    // Declare all vectors/matrices
    real vector alldata, grid_rand, target_cdf, w, sol
    real matrix control_cdf, c, ecmat, lowerbd, upperbd
    real vector bec
    
    // Get dimensions
    n = length(target)
    J = cols(controls)
    
    // Combine all data to find support
    alldata = target
    for (i=1; i<=J; i++) {
        alldata = alldata \ controls[,i]
    }
    
    // Set up random grid
    grid_min = min(alldata)
    grid_max = max(alldata)
    grid_rand = runiform(M,1)*(grid_max - grid_min) + J(M,1, grid_min)
    
    // Compute CDFs at random points
    target_cdf = cdf_at_points(target, grid_rand)
    control_cdf = J(M,J,.)
    for (i=1; i<=J; i++) {
        control_cdf[,i] = cdf_at_points(controls[,i], grid_rand)
    }
	
    
    // Set up LP problem
    // Variables: [w_1,...,w_J,s_1+,...,s_M+,s_1-,...,s_M-]
    total_vars = J + 2*M
    
    // Objective coefficients (minimize sum of s+ and s-)
    c = J(1, total_vars, 0)
    c[1,(J+1)..(J+M)] = J(1,M,1)
    c[1,(J+M+1)..(J+2*M)] = J(1,M,1)
    
    // Equality constraints
    ecmat = J(M+1, total_vars, 0)
    bec = J(M+1,1,0)
    
    // Sum of weights = 1
    ecmat[1,1..J] = J(1,J,1)
    bec[1] = 1
    
    // CDF matching constraints
    for (g=1; g<=M; g++) {
        ecmat[1+g,1..J] = control_cdf[g,.]
        ecmat[1+g,J+g] = -1        // s+
        ecmat[1+g,J+M+g] = 1       // s-
        bec[1+g] = target_cdf[g]
    }
    
    // Set up bounds
    lowerbd = J(1,total_vars,.)
    upperbd = J(1,total_vars,.)
    
    // Weight bounds for simplex case
    if (simplex==1) {
        lowerbd[1,1..J] = J(1,J,0)
    }
    // Non-negative bounds for s+ and s-
    lowerbd[1,(J+1)..(J+M)] = J(1,M,0)
    lowerbd[1,(J+M+1)..(J+2*M)] = J(1,M,0)
    
    // Solve LP
    class LinearProgram scalar q
    q = LinearProgram()
    q.setCoefficients(c)
    q.setMaxOrMin("min")
    q.setEquality(ecmat,bec)
    q.setBounds(lowerbd,upperbd)
    
    val = q.optimize()
    if (q.errorcode()!=0) {
        // Return uniform weights if optimization fails
        w = J(J,1,1/J)
        return(w)
    }
    
    // Extract and normalize weights
    sol = q.parameters()
    w = sol[1..J]'
    ss = sum(w)
    if (abs(ss-1)>1e-8) w = w/ss
    
    return(w)
}


// Compute CDF values at specified grid points
// Parameters:
//   x: Input vector
//   grid: Vector of points at which to evaluate CDF
real vector cdf_builder(real vector x, real vector grid) {
    real vector Xs, out
    real scalar N, G, g, pos
    
    Xs = sort(x, 1)
    N = length(Xs)
    G = length(grid)
    out = J(G,1,0)
    
    // Calculate CDF at each grid point
    for (g=1; g<=G; g++) {
        pos = sum(Xs:<=grid[g])
        out[g] = pos/N
    }
    return(out)
}

// Helper to compute CDF at given grid points
real vector cdf_at_points(real vector x, real vector grid) {
    real scalar N, G, g, pos
    real vector xs, out
    N = length(x)
    xs = sort(x, 1)
    G = length(grid)
    out = J(G, 1, 0)
    for (g=1; g<=G; g++) {
        pos = sum(xs :<= grid[g])
        out[g] = pos / N
    }
    return(out)
}

struct disco_out {
    real matrix weights, quantile_diff, cdf_diff, quantile_synth, cdf_synth,  
	quantile_t, cdf_t, cids
}
 // Main DISCO (DIStribution COmparison) function implementation
// Parameters:
//   y: real vector - Outcome variable for all units and time periods
//   id: real vector - Unit identifiers matching y values
//   tt: real vector - Time period indicators matching y values
//   target_id: real scalar - Identifier for the treated unit
//   T0: real scalar - First post-treatment period
//   T_max: real scalar - Total number of time periods
//   M: real scalar - Number of quantile points
//   G: real scalar - Number of grid points for CDF
//   q_min: real scalar - Minimum quantile level (e.g., 0.1)
//   q_max: real scalar - Maximum quantile level (e.g., 0.9)
//   simplex: real scalar - If 1, weights are non-negative
//   mixture: real scalar - If 1, uses CDF-based method
struct disco_out disco_full_run(real vector y, real vector id, real vector tt,
                              real scalar target_id, real scalar T0, real scalar T_max, 
                              real scalar M, real scalar G, real scalar q_min, real scalar q_max,
                              real scalar simplex, real scalar mixture) 
{
    // Declare all vectors
    real vector uid, cids, yt, idt, target_data, cd, w, w_temp, W_avg, 
             Q_synth, Q_synth_sorted, cids_t, grid, Tq, Tc, C_synth
    
    // Declare all matrices
    real matrix controls_data, Cq, Cc, quantile_diff, cdf_diff, weights_store, 
             period_weights, Q_target_all, Q_synth_all, C_target_all, C_synth_all
    
    // Declare all scalars
    real scalar J_sc, amin, amax, gg, pos, t, ci, m, p, ci2, gg2, gg3, m2, J_t, ci5, ci6, ci7, t2

    // Initialize study units
    uid = unique(id)
    cids = select(uid, uid:!=target_id)
    J_sc = length(cids)

    // Set up evaluation grid
    amin = min(y)
    amax = max(y)
    grid = range(amin, amax, (amax - amin) / (G - 1))'
	
	// Store grid bounds in e()
    st_numscalar("amin", amin)
    st_numscalar("amax", amax)

    // Initialize storage matrices
    quantile_diff = J(M,T_max,.)
    cdf_diff = J(G,T_max,.)
    weights_store = J(T0-1,J_sc,.)
    period_weights = J(T_max,J_sc,.)
    Q_target_all = J(M,T_max,.)
    Q_synth_all = J(M,T_max,.)
    C_target_all = J(G,T_max,.)
    C_synth_all = J(G,T_max,.)

    // First pass: compute weights for pre-treatment period
    for (t=1; t<=T_max; t++) {
        // Extract period data
        yt = select(y, tt:==t)
        idt = select(id, tt:==t)
        target_data = select(yt, idt:==target_id)
        controls_data = .

        // Construct control data matrix
        for (ci2=1; ci2<=J_sc; ci2++) {
            cd = select(yt, idt:==cids[ci2])
            if (ci2==1) {
                controls_data = cd
            }
            else {
                controls_data = controls_data, cd
            }
        }

        // Compute target distributions
        Tq = disco_quantile(target_data,M,q_min,q_max)
        Tc = cdf_builder(target_data,grid)

        // Compute weights for pre-treatment period
        if (t<=T0-1) {
            if (mixture==0) {
                w = disco_solve_weights(controls_data, target_data, M, q_min, q_max, "", 7, simplex)
            } else {
                w_temp = disco_mixture_weights(controls_data, target_data, M, simplex)
                w = w_temp' // transpose
            }
            weights_store[t,.] = w
            period_weights[t,.] = w
        } else {
            period_weights[t,.] = J(J_sc,1,.)'
        }

        Q_target_all[,t] = Tq
        C_target_all[,t] = Tc
    }

    // Compute average weights for post-treatment period
    W_avg = (colsum(weights_store)/((T0-1)))'
    for (t2=T0; t2<=T_max; t2++) {
        period_weights[t2,.] = W_avg'
    }

    // Second pass: compute synthetic controls using average weights
    for (t=1; t<=T_max; t++) {
        // Extract period data
        yt = select(y, tt:==t)
        idt = select(id, tt:==t)
        target_data = select(yt, idt:==target_id)
        cids_t = select(unique(idt), unique(idt):!=target_id)
        
        // Construct control data
        J_t = length(cids_t)
        controls_data = .
        for (ci5=1; ci5<=J_t; ci5++) {
            cd = select(yt, idt:==cids_t[ci5])
            if (ci5==1) controls_data = cd
            else controls_data = controls_data, cd
        }

        w = period_weights[t,.]'

        if (mixture==0) {
            // Quantile-based synthetic control
            real matrix Cq2
            Cq2 = J(M,J_t,.)
            for (ci6=1; ci6<=J_t; ci6++) {
                Cq2[,ci6] = disco_quantile(controls_data[,ci6],M,q_min,q_max)
            }
            Q_synth = Cq2*w
            Q_synth_sorted = sort(Q_synth, 1)
            
            // Calculate synthetic CDF from Q_synth
            C_synth = J(G,1,0)
            for (gg2=1; gg2<=G; gg2++) {
                pos = sum(Q_synth_sorted:<=grid[gg2])
                C_synth[gg2] = pos/M
            }
        } else {
            // CDF-based synthetic control
            real matrix Cc2
            Cc2 = J(G,J_t,.)
            for (ci7=1; ci7<=J_t; ci7++) {
                Cc2[,ci7] = cdf_builder(controls_data[,ci7],grid)
            }
            C_synth = Cc2*w
            Q_synth = J(M,1,.)

            // Calculate synthetic quantile from CDF_synth
            for (m2=1; m2<=M; m2++) {
                p = (m2-1)/(M-1)
                gg3 = 1
                while (gg3<G && C_synth[gg3]<p) {
                    gg3++
                }
                if (gg3>G) gg3=G
                Q_synth[m2] = grid[gg3]
            }
        }
        
        // Store results
        Q_synth_all[,t] = Q_synth
        C_synth_all[,t] = C_synth
        quantile_diff[,t] = Q_target_all[,t] - Q_synth
        cdf_diff[,t] = C_target_all[,t] - C_synth
    }

    // Return results
    struct disco_out scalar r
    r.weights = W_avg
    r.quantile_diff = quantile_diff
    r.cdf_diff = cdf_diff
    r.quantile_synth = Q_synth_all
    r.cdf_synth = C_synth_all
    r.quantile_t = Q_target_all
    r.cdf_t = C_target_all
	r.cids = cids
    
    return(r)
}



// Compute ratio of post/pre treatment effects for placebo testing
// Parameters:
//   y: Outcome variable
//   id: Unit identifiers
//   tt: Time period indicators
//   target_id: ID of treated unit
//   T0: Pre-treatment period
//   T_max: Total number of periods
//   M, G: Number of quantile and grid points
//   q_min, q_max: Range for quantile computation
//   simplex: Whether to constrain weights to be non-negative
//   mixture: Whether to use mixture method
real scalar disco_compute_ratio(real vector y, id, tt, 
                              real scalar target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture) 
{
    // Run main DISCO implementation
    struct disco_out scalar rr
    rr = disco_full_run(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)
    
    // Declare aggregation variables
    real scalar pre_dist, pre_count, post_dist, post_count, dist_t, ratio, t
    
    // Initialize counters
    pre_dist = pre_count = post_dist = post_count = 0
    
    // Calculate mean squared quantile differences for pre/post periods
    for (t=1; t<=T_max; t++) {
        dist_t = mean((rr.quantile_diff[,t]:^2))
        if (t<T0) {
            pre_dist = pre_dist + dist_t
            pre_count = pre_count + 1
        } else {
            post_dist = post_dist + dist_t
            post_count = post_count + 1
        }
    }
    
    // Return missing if insufficient data
    if (pre_count==0 | post_count==0) return(.)
    
    // Compute ratio of post/pre RMSE
    ratio = sqrt((post_dist/post_count))/sqrt((pre_dist/pre_count))
    return(ratio)
}

// Permutation test implementation
// Computes test statistic for each control unit and returns p-value
real scalar disco_permutation_test(real vector y, id, tt,
                                 real scalar target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture) 
{
    // Declare variables
    real scalar actual_ratio, pval, rj, J, count, j
    real vector uid, cids
    
    // Compute ratio for treated unit
    actual_ratio = disco_compute_ratio(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)
    
    // Get control units
    uid = unique(id)
    cids = select(uid, uid:!=target_id)
    J = length(cids)
    
    // Count number of control units with larger ratio
    count = 0
    for (j=1; j<=J; j++) {
        rj = disco_compute_ratio(y, id, tt, cids[j], T0, T_max, M, G, q_min, q_max, simplex, mixture)
        if (rj>=actual_ratio) count = count + 1
    }
    
    // Compute p-value
    pval = (count+1)/(J+1)
    return(pval)
}

// Output structure for confidence intervals
struct CI_out {
    real matrix qdiff_lower,    // Lower bound for quantile differences
              qdiff_upper,      // Upper bound for quantile differences
              cdiff_lower,      // Lower bound for CDF differences
              cdiff_upper       // Upper bound for CDF differences
}

// Output structure for bootstrap iteration
struct iter_out {
    real vector target_q,       // Target quantiles
              target_cdf,       // Target CDF
              weights          // Optimal weights
    real matrix controls_q,     // Control unit quantiles
              controls_cdf      // Control unit CDFs
}

// Output structure for bootstrap results
struct boot_out {
    real matrix quantile_diff,  // Quantile differences
              cdf_diff,        // CDF differences
              quantile_synth,  // Synthetic quantiles
              cdf_synth,      // Synthetic CDFs
              quantile_t,     // Target quantiles
              cdf_t          // Target CDFs
}

/**
 * Computes confidence interval bounds for quantile and CDF differences
 * from bootstrap samples. Can compute either pointwise intervals or uniform bands.
 *
 * Parameters
 * ----------
 * quantile_diff_boot : real matrix
 *     Matrix of bootstrap quantile differences. Each column is one bootstrap sample.
 *     Dimensions: (M*T_max) × boots, where M is number of quantile points,
 *     T_max is number of time periods, boots is number of bootstrap samples
 * 
 * cdf_diff_boot : real matrix
 *     Matrix of bootstrap CDF differences. Each column is one bootstrap sample.
 *     Dimensions: (G*T_max) × boots, where G is number of grid points
 * 
 * M : real scalar
 *     Number of quantile points
 * 
 * G : real scalar
 *     Number of grid points for CDF
 * 
 * T_max : real scalar
 *     Total number of time periods
 * 
 * cl : real scalar
 *     Confidence level (e.g., 0.95 for 95% confidence intervals)
 * 
 * uniform : real scalar
 *     If 1, computes uniform confidence bands
 *     If 0, computes pointwise confidence intervals
 * 
 * main_run : struct disco_out scalar
 *     Results from main DISCO analysis, needed for uniform bands
 * 
 * Returns
 * -------
 * struct CI_out scalar containing:
 *     qdiff_lower : M × T_max matrix of lower bounds for quantile differences
 *     qdiff_upper : M × T_max matrix of upper bounds for quantile differences
 *     cdiff_lower : G × T_max matrix of lower bounds for CDF differences
 *     cdiff_upper : G × T_max matrix of upper bounds for CDF differences
 */
struct CI_out scalar compute_CI_bounds(real matrix quantile_diff_boot, cdf_diff_boot,
                                     real scalar M, G, T_max, cl, uniform,
                                     struct disco_out scalar main_run)
{
    struct CI_out scalar co
    
    // Initialize CI matrices
    co.qdiff_lower = J(M, T_max, .)
    co.qdiff_upper = J(M, T_max, .)
    co.cdiff_lower = J(G, T_max, .)
    co.cdiff_upper = J(G, T_max, .)
    
    if(!uniform) {
        // Compute pointwise confidence intervals
        real scalar alpha, lower_idx, upper_idx, idx
		real vector vals

        
        // Calculate indices for confidence bounds based on confidence level
        alpha = (1-cl)/2  // e.g., 0.025 for 95% CI
        lower_idx = max((1, ceil(alpha*cols(quantile_diff_boot))))
        upper_idx = min((cols(quantile_diff_boot), ceil((1-alpha)*cols(quantile_diff_boot))))
        
        // Process each time period
        real scalar t
        for(t=1; t<=T_max; t++) {
            // Compute quantile difference bounds
            for(idx=1; idx<=M; idx++) {
                // Extract and sort bootstrap values for this quantile and time
                vals = quantile_diff_boot[(t-1)*M + idx,.]'
                vals = sort(vals, 1)
                
                // Get confidence bounds from sorted values
                co.qdiff_lower[idx,t] = vals[lower_idx]
                co.qdiff_upper[idx,t] = vals[upper_idx]
            }
            
            // Compute CDF difference bounds
            for(idx=1; idx<=G; idx++) {
                // Extract and sort bootstrap values for this grid point and time
                vals = cdf_diff_boot[(t-1)*G + idx,.]'
                vals = sort(vals, 1)
                
                // Get confidence bounds from sorted values
                co.cdiff_lower[idx,t] = vals[lower_idx]
                co.cdiff_upper[idx,t] = vals[upper_idx]
            }
        }
    }
    else {
        // Compute uniform confidence bands
        
        // Initialize vectors for maximum absolute deviations
        real vector qmax_abs, cmax_abs
		real matrix qdiff_mat, cdiff_mat, qdiff_err, cdiff_err

        qmax_abs = J(cols(quantile_diff_boot), 1, .)
        cmax_abs = J(cols(quantile_diff_boot), 1, .)
        
        // Calculate maximum absolute deviations for each bootstrap sample
        real scalar b
        for(b=1; b<=cols(quantile_diff_boot); b++) {
            
            // Reshape bootstrap sample to original dimensions
            qdiff_mat = rowshape(quantile_diff_boot[,b], M)
            cdiff_mat = rowshape(cdf_diff_boot[,b], G)
            
            // Compute deviations from main analysis
            qdiff_err = qdiff_mat :- main_run.quantile_diff
            cdiff_err = cdiff_mat :- main_run.cdf_diff
            
            // Store maximum absolute deviations
            qmax_abs[b] = max(abs(vec(qdiff_err)))
            cmax_abs[b] = max(abs(vec(cdiff_err)))
        }
        
        // Compute critical values for uniform bands
        real scalar q_crit, c_crit, m_i, g_i, base_val

        real vector tmp
        
        // Get quantile critical value
        tmp = sort(qmax_abs, 1)
        q_crit = tmp[ceil((1-(1-cl)/2)*cols(quantile_diff_boot))]
        
        // Get CDF critical value
        tmp = sort(cmax_abs, 1)
        c_crit = tmp[ceil((1-(1-cl)/2)*cols(quantile_diff_boot))]
        
        // Calculate uniform bands for each time period
        real scalar t2
        for(t2=1; t2<=T_max; t2++) {
            
            // Compute quantile difference bands
            for(m_i=1; m_i<=M; m_i++) {
                base_val = main_run.quantile_diff[m_i,t2]
                co.qdiff_lower[m_i,t2] = base_val - q_crit
                co.qdiff_upper[m_i,t2] = base_val + q_crit
            }
            
            // Compute CDF difference bands
            for(g_i=1; g_i<=G; g_i++) {
                base_val = main_run.cdf_diff[g_i,t2]
                co.cdiff_lower[g_i,t2] = base_val - c_crit
                co.cdiff_upper[g_i,t2] = base_val + c_crit
            }
        }
    }
    
    return(co)
}

// Performs a single bootstrap iteration for confidence interval computation
// Parameters:
//   y: real vector - Outcome variable for all units and time periods
//   id: real vector - Unit identifiers matching y values
//   tt: real vector - Time period indicators matching y values
//   target_id: real scalar - Identifier for treated unit
//   t: real scalar - Current time period
//   T0: real scalar - First post-treatment period
//   M: real scalar - Number of quantile points
//   G: real scalar - Number of grid points for CDF
//   grid: real vector - Grid points for CDF evaluation
//   q_min: real scalar - Minimum quantile level
//   q_max: real scalar - Maximum quantile level
//   simplex: real scalar - If 1, weights are non-negative
//   mixture: real scalar - If 1, uses CDF-based method
struct iter_out disco_CI_iter(real vector y, real vector id, real vector tt,
                            real scalar target_id, t, T0, M, G,
                            real vector grid,
                            real scalar q_min, q_max, simplex, mixture) 
{
    // Declare structure and variables at start
    struct iter_out scalar out
    real scalar t_len, c_len, J
    real vector yt, idt, target_data, mytar, uid, cids, cd, mycon, indices_t, indices_c
    real matrix mycon_q, mycon_cdf, controls_resampled
    
    // Extract period data once
    yt = select(y, tt:==t)
    idt = select(id, tt:==t)
    target_data = select(yt, idt:==target_id)
    
    // Bootstrap target data
    t_len = length(target_data)
    indices_t = ceil(runiform(t_len,1):*t_len)
    mytar = target_data[indices_t]
    
    // Compute target distributions
    out.target_q = disco_quantile(mytar, M, q_min, q_max)
    out.target_cdf = cdf_builder(mytar, grid)
    
    // Get control units
    uid = unique(idt)
    cids = select(uid, uid:!=target_id)
    J = length(cids)
    
    // Initialize control matrices with known dimensions
    mycon_q = J(M, J, .)
    mycon_cdf = J(G, J, .)
    
    // Bootstrap and compute distributions for each control
    for (ci=1; ci<=J; ci++) {
        // Extract and bootstrap control data
        cd = select(yt, idt:==cids[ci])
        c_len = length(cd)
        indices_c = ceil(runiform(c_len,1):*c_len)
        mycon = cd[indices_c]
        
        // Compute control distributions
        mycon_q[,ci] = disco_quantile(mycon, M, q_min, q_max)
        mycon_cdf[,ci] = cdf_builder(mycon, grid)
    }
    
    // Initialize weights vector
    out.weights = J(J, 1, .)
    
    // Compute weights for pre-treatment period
    if (t<=T0) {
        if (mixture==0) {
            // Quantile-based weights
            // First attempt with single bootstrap sample
            out.weights = disco_solve_weights(mycon, mytar, M, q_min, q_max, "", 7, simplex)
            
            // Assemble matrix with all resampled controls
            controls_resampled = mycon
            if (J>1) {
                for (cc=2; cc<=J; cc++) {
                    cd = select(yt, idt:==cids[cc])
                    c_len = length(cd)
                    indices_c = ceil(runiform(c_len,1):*c_len)
                    mycon = cd[indices_c]
                    controls_resampled = controls_resampled, mycon
                }
            }
            // Recompute weights with full resampled matrix
            out.weights = disco_solve_weights(controls_resampled, mytar, M, q_min, q_max, "", 7, simplex)
        } 
        else {
            // CDF-based weights
            out.weights = disco_mixture_weights(mycon_cdf, out.target_cdf, simplex)
        }
    }
    
    // Store control distributions
    out.controls_q = mycon_q
    out.controls_cdf = mycon_cdf
    
    return(out)
}
// Generate bootstrap counterfactuals for confidence interval computation
// Parameters:
//   iter_results: Vector of iteration results from bootstrap samples
//   T0: Pre-treatment period
//   T_max: Total periods
//   M, G: Number of quantile and grid points
//   grid: Vector of grid points for CDF
//   mixture: Whether to use mixture method
// Update the bootCounterfactuals function
struct boot_out bootCounterfactuals(struct iter_out vector iter_results,
                                  real scalar T0, T_max, M, G,
                                  real vector grid,
                                  real scalar mixture) 
{
    // Declare output structure and variables
    struct boot_out scalar bo
    real scalar J, t, gg, m2, gg3, p
    real vector W_avg, target_q, target_cdf, Q_synth, C_synth, Q_synth_sorted
    real matrix weights_all, quantile_diff, cdf_diff, mycon_q, mycon_cdf
    
    // Get number of control units
    J = cols(iter_results[1].controls_q)
    
    // Initialize storage matrices
    quantile_diff = J(M,T_max,.)
    cdf_diff = J(G,T_max,.)
    quantile_synth = J(M,T_max,.)
    cdf_synth = J(G,T_max,.)
    quantile_t = J(M,T_max,.)
    cdf_t = J(G,T_max,.)
    
    // Compute average weights from pre-treatment period
    weights_all = J(T0-1,J,.)
    for (t=1; t<=T0-1; t++) {
        weights_all[t,.] = (iter_results[t].weights)
    }
    W_avg = (colsum(weights_all)/(T0-1))'
    
    // Compute differences for each time period
    for (t=1; t<=T_max; t++) {
        // Extract target values
        target_q = iter_results[t].target_q
        target_cdf = iter_results[t].target_cdf
        mycon_q = iter_results[t].controls_q
        mycon_cdf = iter_results[t].controls_cdf
        
        // Store target values
        quantile_t[,t] = target_q
        cdf_t[,t] = target_cdf
        
        // Compute synthetic control based on method
        if (mixture==0) {
            // Quantile-based synthetic control
            Q_synth = mycon_q*W_avg
            Q_synth_sorted = sort(Q_synth,1)
            
            // Compute corresponding CDF
            C_synth = J(G,1,0)
            for (gg=1; gg<=G; gg++) {
                C_synth[gg] = sum(Q_synth_sorted:<=grid[gg])/M
            }
        } else {
            // CDF-based synthetic control
            C_synth = mycon_cdf*W_avg
            
            // Compute corresponding quantiles
            Q_synth = J(M,1,.)
            for (m2=1; m2<=M; m2++) {
                p = (m2-1)/(M-1)
                gg3 = 1
                while (gg3<=G & C_synth[gg3]<p) gg3++
                Q_synth[m2] = grid[gg3>G ? G : gg3]
            }
        }
        
        // Store results
        quantile_synth[,t] = Q_synth
        cdf_synth[,t] = C_synth
        quantile_diff[,t] = target_q - Q_synth
        cdf_diff[,t] = target_cdf - C_synth
    }
    
    // Set output fields
    bo.quantile_diff = quantile_diff
    bo.cdf_diff = cdf_diff
    bo.quantile_synth = quantile_synth
    bo.cdf_synth = cdf_synth
    bo.quantile_t = quantile_t
    bo.cdf_t = cdf_t
    return(bo)
}


// Computes bootstrap confidence intervals for DISCO estimates
// Parameters:
//   y: real vector - Outcome variable for all units and time periods
//   id: real vector - Unit identifiers matching y values
//   tt: real vector - Time period indicators matching y values
//   target_id: real scalar - Identifier for treated unit
//   T0: real scalar - First post-treatment period
//   T_max: real scalar - Total number of time periods
//   M: real scalar - Number of quantile points
//   G: real scalar - Number of grid points for CDF
//   q_min, q_max: real scalar - Quantile range bounds
//   simplex: real scalar - If 1, weights are non-negative
//   mixture: real scalar - If 1, uses CDF-based method
//   boots: real scalar - Number of bootstrap iterations
//   cl: real scalar - Confidence level (e.g., 0.95)
//   uniform: real scalar - If 1, compute uniform confidence bands
// Complete bootstrap confidence interval computation function
struct CI_out scalar disco_bootstrap_CI(real vector y, real vector id, real vector tt,
                                      real scalar target_id, T0, T_max, M, G,
                                      real scalar q_min, q_max, simplex, mixture, boots, cl, uniform) 
{
    // Declare structures and core variables
    struct CI_out scalar co
    struct boot_out scalar bo
    struct disco_out scalar main_run
    struct iter_out vector iter_results
    real scalar amin, amax, N, b, t, alpha, lower_idx, upper_idx, idx, b2, q_crit, c_crit
    real vector grid, vals, tmp, qmax_abs, cmax_abs
    real matrix quantile_diff_boot, cdf_diff_boot, qdiff_lower, qdiff_upper, cdiff_lower, cdiff_upper
    
    // Additional matrices for storing bootstrap results
    real matrix quantile_synth_boot, cdf_synth_boot, quantile_t_boot, cdf_t_boot
    
    // Initialize iteration results vector
    iter_results = iter_out(T_max)
    
    // Set up evaluation grid
    amin = min(y)
    amax = max(y)
    grid = range(amin, amax, (amax - amin)/(G-1))'
    
    // Run main analysis
    main_run = disco_full_run(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)
    
    // Initialize bootstrap storage
    N = length(y)
    quantile_diff_boot = J(M*T_max, boots, .)
    cdf_diff_boot = J(G*T_max, boots, .)
    quantile_synth_boot = J(M*T_max, boots, .)
    cdf_synth_boot = J(G*T_max, boots, .)
    quantile_t_boot = J(M*T_max, boots, .)
    cdf_t_boot = J(G*T_max, boots, .)
    
    // Perform bootstrap iterations
    for (b=1; b<=boots; b++) {
        // Compute bootstrap sample for each time period
        for (t=1; t<=T_max; t++) {
            struct iter_out scalar out
            out = disco_CI_iter(y, id, tt, target_id, t, T0, M, G, grid, q_min, q_max, simplex, mixture)
            iter_results[t] = out
        }
        
        // Compute counterfactuals for this bootstrap sample
        bo = bootCounterfactuals(iter_results, T0, T_max, M, G, grid, mixture)
        
        // Store all bootstrap results
        quantile_diff_boot[,b] = vec(bo.quantile_diff)
        cdf_diff_boot[,b] = vec(bo.cdf_diff)
        quantile_synth_boot[,b] = vec(bo.quantile_synth)
        cdf_synth_boot[,b] = vec(bo.cdf_synth)
        quantile_t_boot[,b] = vec(bo.quantile_t)
        cdf_t_boot[,b] = vec(bo.cdf_t)
    }
    
    // Initialize confidence interval matrices
    qdiff_lower = J(M,T_max,.)
    qdiff_upper = J(M,T_max,.)
    cdiff_lower = J(G,T_max,.)
    cdiff_upper = J(G,T_max,.)
    
    if (uniform==0) {
        // Compute pointwise confidence intervals
        for (t=1; t<=T_max; t++) {
            // Process quantile differences
            for (idx=1; idx<=M; idx++) {
                real scalar row_i
                row_i = (t-1)*M + idx
                vals = quantile_diff_boot[row_i,.]'
                tmp = sort(vals, 1)
                
                // Calculate confidence bounds
                alpha = (1-cl)/2
                lower_idx = max((ceil(alpha*boots), 1))
                upper_idx = min((ceil((1-alpha)*boots), boots))
                
                qdiff_lower[idx,t] = tmp[lower_idx]
                qdiff_upper[idx,t] = tmp[upper_idx]
            }
            
            // Process CDF differences
            for (idx=1; idx<=G; idx++) {
                real scalar row_g
                row_g = (t-1)*G + idx
                vals = cdf_diff_boot[row_g,.]'
                tmp = sort(vals, 1)
                
                alpha = (1-cl)/2
                lower_idx = max((ceil(alpha*boots), 1))
                upper_idx = min((ceil((1-alpha)*boots), boots))
                
                cdiff_lower[idx,t] = tmp[lower_idx]
                cdiff_upper[idx,t] = tmp[upper_idx]
            }
        }
    } 
    else {
        // Compute uniform confidence bands
        qmax_abs = J(boots,1,.)
        cmax_abs = J(boots,1,.)
        
        // Calculate maximum absolute deviations
        for (b2=1; b2<=boots; b2++) {
            real matrix qdiff_mat, cdiff_mat, qdiff_err, cdiff_err
            
            qdiff_mat = rowshape(quantile_diff_boot[,b2], M)
            cdiff_mat = rowshape(cdf_diff_boot[,b2], G)
            
            qdiff_err = qdiff_mat - main_run.quantile_diff
            cdiff_err = cdiff_mat - main_run.cdf_diff
            
            qmax_abs[b2] = max(abs(vec(qdiff_err)))
            cmax_abs[b2] = max(abs(vec(cdiff_err)))
        }
        
        // Compute critical values
        tmp = sort(qmax_abs, 1)
        q_crit = tmp[ceil((1-(1-cl)/2)*boots)]
        tmp = sort(cmax_abs, 1)
        c_crit = tmp[ceil((1-(1-cl)/2)*boots)]
        
        // Calculate uniform bands
        for (t=1; t<=T_max; t++) {
            real scalar m_i, g_i, base_val
            for (m_i=1; m_i<=M; m_i++) {
                base_val = main_run.quantile_diff[m_i,t]
                qdiff_lower[m_i,t] = base_val - q_crit
                qdiff_upper[m_i,t] = base_val + q_crit
            }
            for (g_i=1; g_i<=G; g_i++) {
                base_val = main_run.cdf_diff[g_i,t]
                cdiff_lower[g_i,t] = base_val - c_crit
                cdiff_upper[g_i,t] = base_val + c_crit
            }
        }
    }
    
    // Store results
    co.qdiff_lower = qdiff_lower
    co.qdiff_upper = qdiff_upper
    co.cdiff_lower = cdiff_lower
    co.cdiff_upper = cdiff_upper
    
    return(co)
}
// Wrapper function for main DISCO implementation
// Handles Stata integration and matrix storage
real scalar disco_wrapper(real vector y, id, tt,
                         real scalar target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture) 
{
    // Run main DISCO analysis
    struct disco_out scalar results
    results = disco_full_run(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)
    
    // Store results in Stata matrices
    st_matrix("weights", results.weights')
    st_matrix("quantile_diff", results.quantile_diff)
    st_matrix("cdf_diff", results.cdf_diff)
    st_matrix("quantile_synth", results.quantile_synth)
    st_matrix("quantile_t", results.quantile_t)
    st_matrix("cdf_synth", results.cdf_synth)
    st_matrix("cdf_t", results.cdf_t)
	st_matrix("cids", results.cids')
    
    return(0)
}

// Wrapper function for permutation test
// Handles Stata integration
real scalar disco_permutation_wrapper(string scalar y_name, id_name, tt_name,
                                    real scalar target_id, T0, T_max, M, G,
                                    real scalar q_min, q_max, simplex, mixture)
{
    // Note: Implementation would load data from Stata
    // and call disco_permutation_test
    return(0)
}

// Wrapper function for confidence interval computation
real scalar disco_ci_wrapper(real vector y, id, tt,
                           real scalar target_id, T0, T_max, M, G,
                           real scalar q_min, q_max, simplex, mixture, boots, cl, uniform)
{
    // Compute confidence intervals
    struct CI_out scalar results
    results = disco_bootstrap_CI(y, id, tt, target_id, T0, T_max, M, G,
                               q_min, q_max, simplex, mixture,
                               boots, cl, uniform)
    
    // Store difference confidence intervals
    st_matrix("qdiff_lower", results.qdiff_lower)
    st_matrix("qdiff_upper", results.qdiff_upper)
    st_matrix("cdiff_lower", results.cdiff_lower)
    st_matrix("cdiff_upper", results.cdiff_upper)
    
    
    return(0)
}




/*
Helper function to identify indices in grid corresponding to quantile range
Inputs:
    q_start: starting quantile value
    q_end: ending quantile value
    grid_length: length of evaluation grid
Returns:
    2x1 vector containing [start_idx, end_idx]
*/
real vector find_grid_indices(real scalar q_start, real scalar q_end, 
                            real scalar grid_length) 
{
    real scalar idx_start, idx_end
    
    // Ensure values are within bounds
    if (q_start < 0 | q_end > 1) {
        _error("Quantile values must be between 0 and 1")
    }
    
    // Convert quantiles to indices
    idx_start = max((1, ceil(q_start * grid_length)))
    idx_end = min((grid_length, floor(q_end * grid_length)))
    
    return((idx_start \ idx_end))
}

/*
Helper function to compute mean effect over a range
Inputs:
    data: matrix containing effects
    idx_start: starting index
    idx_end: ending index
    t: time period column
Returns:
    Mean effect over specified range
*/
real scalar compute_mean_effect(real matrix data, real scalar idx_start, 
                              real scalar idx_end, real scalar t) 
{
    real vector range
    range = data[|idx_start, t \ idx_end, t|]
    
    // Check for missing values
    if (anyof(range, .)) {
        return(.)
    }
    
    return(mean(range))
}

/*
Helper function to compute bootstrap statistics
Inputs:
    boot_data: 3D array of bootstrap samples
    idx_start: starting index
    idx_end: ending index
    t: time period
    cl: confidence level
Returns:
    3x1 vector containing [se, ci_lower, ci_upper]
*/
real vector compute_boot_stats(real matrix boot_data, real scalar idx_start, 
                             real scalar idx_end, real scalar t, 
                             real scalar cl) 
{
    real vector boot_means, stats
    
    // Compute mean effect for each bootstrap sample
    boot_means = J(cols(boot_data), 1, .)
    for(b=1; b<=cols(boot_data); b++) {
        boot_means[b] = mean(boot_data[|idx_start, t, b \ idx_end, t, b|])
    }
    
    // Remove missing values
    boot_means = select(boot_means, !missing(boot_means))
    
    if (length(boot_means) == 0) {
        return(J(3, 1, .))
    }
    
    // Compute standard error
    real scalar se
    se = sqrt(variance(boot_means))
    
    // Compute confidence interval bounds
    real vector ci_bounds
    ci_bounds = get_ci_bounds(boot_means, cl)
    
    return((se \ ci_bounds))
}



/*
Main function to compute summary statistics
Inputs:
    agg: string indicating type of aggregation ("quantile", "cdf", "quantileDiff", "cdfDiff")
    sample_points: vector of quantile points for summary statistics
    T0: first treatment period
    T_max: last period
	quantile_diff: estimates of quantile effects
	cdf_diff: estimates of cdf effects
	CI: 1-0 whether or not to produce confidence intervals
Returns:
    0 if successful, error code otherwise
*/
real scalar compute_summary_stats(string scalar agg, real vector sample_points, 
                                real scalar T0, real scalar T_max, real matrix quantile_diff,
                                real matrix cdf_diff, real scalar CI) 
{
    // Input validation
    if (!anyof(("quantile", "cdf", "quantileDiff", "cdfDiff"), agg)) {
        errprintf("Invalid aggregation type\n")
        return(1)
    }
    
    // Determine if we're doing CDF-based analysis
    real scalar is_cdf
    is_cdf = (agg == "cdf" | agg == "cdfDiff")
    
    // Process grid points based on type
    real vector grid_points
    if (is_cdf) {
        // Get Y support bounds
        real scalar amin, amax
        amin = st_numscalar("e(amin)")
        amax = st_numscalar("e(amax)")
        
        // Convert probability points to Y values if needed
        if (max(sample_points) <= 1 & min(sample_points) >= 0) {
            grid_points_temp = amin :+ sample_points :* (amax - amin)
            // Add endpoints if not included
			amin 
			amax
            if (min(grid_points_temp) > amin) grid_points = amin \ grid_points_temp'
            if (max(grid_points_temp) < amax) grid_points = grid_points \ amax
        } else {
            grid_points = sort(sample_points', 1)
        }
    } else {
        // For quantiles, use probability points directly
        grid_points = sort(sample_points', 1)
        // Add endpoints if not included
        if (min(grid_points) > 0) grid_points = 0 \ grid_points
        if (max(grid_points) < 1) grid_points = grid_points \ 1
    }
    
    // Number of intervals
    real scalar n_intervals

    n_intervals = length(grid_points) - 1
    
    // Initialize summary stats matrix
    real matrix summary_stats
    summary_stats = J(n_intervals * (T_max - T0 + 1), 7, .)
    
    // Row counter
    real scalar row, G, M
    row = 1
	real vector grid, idx, prob_grid, boot_stats
    
    // For each post-treatment period
    real scalar t, i
    for(t = T0; t <= T_max; t++) {
        // For each interval
        for(i = 1; i <= n_intervals; i++) {
            // Store interval bounds and time
            summary_stats[row, 1] = t
            summary_stats[row, 2] = grid_points[i]
            summary_stats[row, 3] = grid_points[i + 1]
            
            // Find indices for this interval
            if (is_cdf) {
                // For CDF, use actual Y values to find relevant grid points
                G = rows(cdf_diff)
                amin = st_numscalar("e(amin)")
                amax = st_numscalar("e(amax)")
                grid = range(amin, amax, (amax - amin)/(G-1))'
                idx = selectindex(grid :>= grid_points[i] :& grid :<= grid_points[i + 1])
            } else {
                // For quantiles, use probability points
                M = rows(quantile_diff)
                prob_grid = range(0, 1, 1/(M-1))'
                idx = selectindex(prob_grid :>= grid_points[i] :& prob_grid :<= grid_points[i + 1])
            }
            
            // Compute effect
            if (length(idx) > 0) {
                if (is_cdf) {
                    summary_stats[row, 4] = mean(cdf_diff[idx, t])
                } else {
                    summary_stats[row, 4] = mean(quantile_diff[idx, t])
                }
            }
            
            // Add confidence intervals if requested
            if (CI) {
                boot_stats = compute_boot_stats(data_boot, idx[1], idx[2], t, cl)
                summary_stats[row, 5] = boot_stats[1]  // SE
                summary_stats[row, 6] = boot_stats[2]  // CI lower
                summary_stats[row, 7] = boot_stats[3]  // CI upper
            }
            
            row++
        }
    }
    
    // Store results in Stata
    st_matrix("summary_stats", summary_stats)
    
    return(0)
}


/*
Helper function to compute confidence interval bounds
Inputs:
    boot_samples: vector of bootstrap samples
    cl: confidence level
Returns:
    2x1 vector containing [lower_bound, upper_bound]
*/
real vector get_ci_bounds(real vector boot_samples, real scalar cl) 
{
    real vector sorted, bounds
    real scalar n, lower_idx, upper_idx
    
    // Sort bootstrap samples
    sorted = sort(boot_samples, 1)
    n = length(sorted)
    
    // Calculate indices for confidence bounds
    lower_idx = ceil((1-cl)/2 * n)
    upper_idx = ceil((1-(1-cl)/2) * n)
    
    // Handle edge cases
    if (lower_idx < 1) lower_idx = 1
    if (upper_idx > n) upper_idx = n
    
    bounds = sorted[lower_idx] \ sorted[upper_idx]
    
    return(bounds)
}

end
