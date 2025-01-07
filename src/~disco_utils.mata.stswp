mata:

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
              cdf_diff,         // CDF differences
              quantile_synth,   // Synthetic quantiles
              cdf_synth,        // Synthetic CDFs
              quantile_t,       // Target quantiles
              cdf_t             // Target CDFs
}

struct disco_out {
    real matrix weights, quantile_diff, cdf_diff, quantile_synth, cdf_synth,  
    quantile_t, cdf_t, cids
}

// Returns get_unique values from a vector
real vector get_unique(real vector x) {
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

    real scalar i
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
real vector disco_quantile(real vector X, real scalar G, real scalar q_min, real scalar q_max) {
    real vector p
    real scalar j, p_j
    
    p = J(G,1,.)
    real scalar j2
    for (j2=1; j2<=G; j2++) {
        p_j = q_min + (q_max - q_min)*(j2-1)/(G-1)
        p[j2] = p_j
    }
    
    return(disco_quantile_points(X,p))
}



// Compute CDF values at specified grid points
real vector cdf_builder(real vector x, real vector grid) {
    real vector Xs, out
    real scalar N, G, g, pos

    Xs = sort(x, 1)
    N = length(Xs)
    G = length(grid)
    out = J(G,1,0)

    real scalar g3
    for (g3=1; g3<=G; g3++) {
        pos = sum(Xs:<=grid[g3])
        out[g3] = pos/N
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

    real scalar g4
    for (g4=1; g4<=G; g4++) {
        pos = sum(xs :<= grid[g4])
        out[g4] = pos / N
    }
    return(out)
}

// ------------------------------------------------------------
// 1. Quadratic-Programming Weights from PRE-COMPUTED Quantiles
// ------------------------------------------------------------
real vector disco_solve_weights(real matrix control_quantiles, real vector target_quantiles, real scalar simplex)
{
    real scalar J
    real matrix Gm, CE, CI, C_in
    real vector g0, ce0, ci0, w, res

    J     = cols(control_quantiles)
    C_in  = control_quantiles

    // Original algebra from disco_solve_weights
    Gm = 2 :* (C_in' * C_in)
    g0 = -2 :* (C_in' * target_quantiles)

    // Sum of weights = 1
    CE  = J(J,1,1)
    ce0 = -1

    // Nonnegative constraints?
    if (simplex==1) {
        CI  = I(J)
        ci0 = J(J,1,0)
    }
    else {
        CI  = J(J,1,1)
        ci0 = -1e20
    }

    // Solve using your existing solver
    res = solve_quadprog(Gm, g0, CE, ce0, CI, ci0)
    w   = res[1..J]

    return(w)
}

// ------------------------------------------------------
// 2. Linear-Programming Weights from PRE-COMPUTED CDFs
// ------------------------------------------------------
real vector disco_mixture_weights(
        real matrix control_cdf,  // G x J matrix of each control's CDF
        real vector target_cdf,   // G x 1 vector of target's CDF
        real scalar simplex       // 1 => simplex constraints, 0 => no constraints
    )
{
    real scalar J, Grows, total_vars
    real matrix c, ecmat, lowerbd, upperbd
    real vector bec, sol, w

    // "Grows" = number of grid points, same as rows in control_cdf
    // J       = number of controls (columns)
    Grows = rows(control_cdf)
    J     = cols(control_cdf)
    total_vars = J + 2*Grows

    //-----------------------------------------------------
    // c is 1×(J + 2*Grows). We need to fill two sub-slices 
    // each of length Grows with 1s. Must assign a 1×Grows matrix.
    //-----------------------------------------------------
    c = J(1, total_vars, 0)
    c[1, (J+1)..(J+Grows)]          = J(1, Grows, 1)    // was "= 1"
    c[1, (J+Grows+1)..(J+2*Grows)]  = J(1, Grows, 1)    // was "= 1"

    //-----------------------------------------------------
    // ecmat is (Grows+1)×(J + 2*Grows). We want the first row 
    // from col 1..J to be 1. Must assign 1×J matrix
    //-----------------------------------------------------
    ecmat = J(Grows+1, total_vars, 0)
    bec   = J(Grows+1, 1, 0)

    // sum of weights = 1
    ecmat[1, 1..J] = J(1, J, 1)  // was "= 1"
    bec[1]         = 1

    //-----------------------------------------------------
    // Build constraints for each grid point 
    //-----------------------------------------------------
    real scalar g2
    for (g2=1; g2<=Grows; g2++) {
        // ecmat[1+g2, 1..J] is a (1×J) slice = control_cdf[g2, .] (1×J)
        ecmat[1+g2, 1..J]       = control_cdf[g2, .]
        // single cells are fine with scalar assignment
        ecmat[1+g2, J+g2]       = -1
        ecmat[1+g2, J+Grows+g2] =  1

        // target_cdf[g2] is also fine as a scalar
        bec[1+g2] = target_cdf[g2]
    }

    //-----------------------------------------------------
    // Bounds
    //-----------------------------------------------------
    lowerbd = J(1, total_vars, .)
    upperbd = J(1, total_vars, .)

    if (simplex == 1) {
        // Nonnegative weights
        lowerbd[1, 1..J] = J(1, J, 0)  // was "= 0"
    }

    // The slacks must be >=0
    lowerbd[1, (J+1)..(J+Grows)]         = J(1, Grows, 0)
    lowerbd[1, (J+Grows+1)..(J+2*Grows)] = J(1, Grows, 0)

    //-----------------------------------------------------
    // Solve LP
    //-----------------------------------------------------
    class LinearProgram scalar q
    q = LinearProgram()
    q.setCoefficients(c)
    q.setMaxOrMin("min")
    q.setEquality(ecmat, bec)
    q.setBounds(lowerbd, upperbd)

    real scalar val
    val = q.optimize()
    if (q.errorcode() != 0) {
        // fallback uniform
        w = J(J, 1, 1/J)
        return(w)
    }
    sol = q.parameters()
    w   = sol[1..J]'   // first J parameters = weights

    //-----------------------------------------------------
    // Enforce sum(w)=1 if there's minor numeric drift
    //-----------------------------------------------------
    real scalar ss
    ss = sum(w)
    if (abs(ss - 1) > 1e-8) w = w / ss

    return(w)
}


// Main DISCO function
struct disco_out disco_full_run(real vector y, real vector id, real vector tt,
                              real scalar target_id, real scalar T0, real scalar T_max, 
                              real scalar G, real scalar q_min, real scalar q_max,
                              real scalar simplex, real scalar mixture) 
{
    real vector uid, cids, yt, idt, target_data, cd, w, w_temp, W_avg, 
             Q_synth, Q_synth_sorted, cids_t, grid, Tq, Tc, C_synth

    // We will create two big storage matrices:
    //   storeCq : (G*T_max) x J_sc  for control quantiles
    //   storeCc : (G*T_max) x J_sc  for control CDF
    real matrix storeCq, storeCc

    real matrix controls_data, Cq, Cc, quantile_diff, cdf_diff, weights_store, 
             period_weights, Q_target_all, Q_synth_all, C_target_all, C_synth_all

    real scalar J_sc, amin, amax, t, ci2, t_loop

    //-------------------------------------------------------
    // 1. Setup: Unique control IDs, define grid, allocate 
    //-------------------------------------------------------
    uid  = get_unique(id)
    cids = select(uid, uid:!=target_id)
    J_sc = length(cids)

    amin = min(y)
    amax = max(y)
    grid = range(amin, amax, (amax - amin) / (G - 1))'

    st_numscalar("amin", amin)
    st_numscalar("amax", amax)

    // For final output
    quantile_diff = J(G, T_max, .)
    cdf_diff      = J(G, T_max, .)
    weights_store = J(T0-1, J_sc, .)
    period_weights = J(T_max, J_sc, .)
    Q_target_all  = J(G, T_max, .)
    Q_synth_all   = J(G, T_max, .)
    C_target_all  = J(G, T_max, .)
    C_synth_all   = J(G, T_max, .)

    // For storing each period's control quantiles/CDF
    // We'll store them rowwise: block of G rows per period
    storeCq = J(G * T_max, J_sc, .)
    storeCc = J(G * T_max, J_sc, .)


    //-------------------------------------------------------
    // 2. First loop: compute quantiles/CDF & weights
    //-------------------------------------------------------
    for (t_loop=1; t_loop<=T_max; t_loop++) {
        // Pull out the data for this period
        yt  = select(y,  tt:==t_loop)
        idt = select(id, tt:==t_loop)

        // Target data
        target_data = select(yt, idt:==target_id)

        // 1. Compute target quantiles, target CDF
        Tq = disco_quantile(target_data, G, q_min, q_max)
        Tc = cdf_builder(target_data, grid)

        // 2. Compute each control's quantiles, control's CDF
        Cq = J(G, J_sc, .)
        Cc = J(G, J_sc, .)

        real scalar ci2_loop
        for (ci2_loop=1; ci2_loop<=J_sc; ci2_loop++) {
            cd = select(yt, idt:==cids[ci2_loop])
            // each column is that control's quantile or CDF
            Cq[,ci2_loop] = disco_quantile(cd, G, q_min, q_max)
            Cc[,ci2_loop] = cdf_builder(cd, grid)
        }

        // 3. Store them in the big containers
        //    So the block for period t_loop is from (t_loop-1)*G+1 to t_loop*G
        storeCq[(t_loop-1)*G+1 .. t_loop*G,  ] = Cq
        storeCc[(t_loop-1)*G+1 .. t_loop*G,  ] = Cc

        // 4. If pre-treatment, solve for weights
        if (t_loop <= T0-1) {
            if (mixture == 0) {
                // Solve from precomputed quantiles
                w = disco_solve_weights(Cq, Tq, simplex)
            } else {
                // Solve from precomputed CDF
                w_temp = disco_mixture_weights(Cc, Tc, simplex)
				w = w_temp'
            }

            weights_store[t_loop, .]  = w
            period_weights[t_loop, .] = w

        }
        else {
            // no weighting in post period
            period_weights[t_loop, .] = J(J_sc,1,.)'
        }

        // 5. Save the target unit's Q & C in Q_target_all, C_target_all
        Q_target_all[, t_loop] = Tq
        C_target_all[, t_loop] = Tc
    }

    //-------------------------------------------------------
    // 3. Compute average weights for post period
    //-------------------------------------------------------
    real scalar T0_minus_1
    T0_minus_1 = T0 - 1
    // If T0-1=0, we should be careful not to divide by zero
    // but presumably T0>1 in your setup
    W_avg = (colsum(weights_store) / T0_minus_1)'

    // Fill in post-period weights with W_avg
    real scalar t2_loop
    for (t2_loop=T0; t2_loop<=T_max; t2_loop++) {
        period_weights[t2_loop, .] = W_avg'
    }

    //-------------------------------------------------------
    // 4. Second loop: use stored Cq/Cc + final weights
    //-------------------------------------------------------
    for (t_loop2=1; t_loop2<=T_max; t_loop2++) {

        // The target Q & C are already in Q_target_all, C_target_all
        // for period t_loop2
        Tq = Q_target_all[, t_loop2]
        Tc = C_target_all[, t_loop2]

        // Retrieve the block of G rows for this period
        real matrix Cq2, Cc2
        Cq2 = storeCq[(t_loop2-1)*G+1 .. t_loop2*G,  ]
        Cc2 = storeCc[(t_loop2-1)*G+1 .. t_loop2*G,  ]

        // get the final weights
        w = period_weights[t_loop2, .]'

        // Build synthetic distribution from stored Q or C
        if (mixture == 0) {
            // Synthesize quantiles
            Q_synth = Cq2 * w
            // Then sort to get synthetic CDF
            Q_synth_sorted = sort(Q_synth, 1)

            C_synth = J(G,1,0)
            real scalar gg2_loop, pos
            for (gg2_loop=1; gg2_loop<=G; gg2_loop++) {
                pos = sum(Q_synth_sorted :<= grid[gg2_loop])
                C_synth[gg2_loop] = pos / G
            }

        } else {
            // Synthesize CDF

            C_synth = Cc2 * w
            // Then invert to get synthetic quantiles
            Q_synth = J(G,1,.)
            real scalar m2_loop, p, gg3
            for (m2_loop=1; m2_loop<=G; m2_loop++) {
                p   = (m2_loop-1)/(G-1)
                gg3 = 1
                while (gg3 < G && C_synth[gg3] < p) {
                    gg3++
                }
                if (gg3 > G) gg3 = G
                Q_synth[m2_loop] = grid[gg3]
            }
        }

        // Fill in the final results for period t_loop2
        Q_synth_all[, t_loop2] = Q_synth
        C_synth_all[, t_loop2] = C_synth

        // Differences
        quantile_diff[, t_loop2] = Tq - Q_synth
        cdf_diff[, t_loop2]      = Tc - C_synth
    }

    //-------------------------------------------------------
    // 5. Return results
    //-------------------------------------------------------
    struct disco_out scalar r
    r.weights         = W_avg
    r.quantile_diff   = quantile_diff
    r.cdf_diff        = cdf_diff
    r.quantile_synth  = Q_synth_all
    r.cdf_synth       = C_synth_all
    r.quantile_t      = Q_target_all
    r.cdf_t           = C_target_all
    r.cids            = cids

    return(r)
}




// Compute ratio
real scalar disco_compute_ratio(real vector y, id, tt, 
                              real scalar target_id, T0, T_max, G, q_min, q_max, simplex, mixture) 
{
    struct disco_out scalar rr
    rr = disco_full_run(y, id, tt, target_id, T0, T_max, G, q_min, q_max, simplex, mixture)
    
    real scalar pre_dist, pre_count, post_dist, post_count, dist_t, ratio, t
    pre_dist = pre_count = post_dist = post_count = 0

    real scalar t_loop
    for (t_loop=1; t_loop<=T_max; t_loop++) {
        dist_t = mean((rr.quantile_diff[,t_loop]:^2))
        if (t_loop<T0) {
            pre_dist = pre_dist + dist_t
            pre_count = pre_count + 1
        } else {
            post_dist = post_dist + dist_t
            post_count = post_count + 1
        }
    }

    if (pre_count==0 | post_count==0) return(.)

    ratio = sqrt((post_dist/post_count))/sqrt((pre_dist/pre_count))
    return(ratio)
}

// Permutation test
real scalar disco_permutation_test(real vector y, id, tt,
                                 real scalar target_id, T0, T_max, G, q_min, q_max, simplex, mixture) 
{
    real scalar actual_ratio, pval, rj, J, count, j
    actual_ratio = disco_compute_ratio(y, id, tt, target_id, T0, T_max, G, q_min, q_max, simplex, mixture)
    
    real vector uid, cids
    uid = get_unique(id)
    cids = select(uid, uid:!=target_id)
    J = length(cids)

    count = 0
    real scalar j_loop
    for (j_loop=1; j_loop<=J; j_loop++) {
        rj = disco_compute_ratio(y, id, tt, cids[j_loop], T0, T_max, G, q_min, q_max, simplex, mixture)
        if (rj>=actual_ratio) count = count + 1
    }

    pval = (count+1)/(J+1)
    return(pval)
}

struct iter_out disco_CI_iter(real vector y, real vector id, real vector tt,
                            real scalar target_id, t, T0, G,
                            real vector grid,
                            real scalar q_min, q_max, simplex, mixture) 
{
    struct iter_out scalar out
    real scalar t_len, c_len, J
    real scalar ci
    real vector yt, idt, target_data, mytar, uid, cids, cd, mycon, indices_t, indices_c
    real matrix mycon_q, mycon_cdf

    //------------------------------------------------------------------
    // 1. Draw a bootstrap sample for the *target* unit in period t
    //------------------------------------------------------------------
    yt = select(y, tt:==t)
    idt = select(id, tt:==t)
    target_data = select(yt, idt:==target_id)

    t_len = length(target_data)
    indices_t = ceil(runiform(t_len,1):*t_len)
    mytar = target_data[indices_t]

    // Compute the target's quantiles & CDF from that bootstrap sample
    out.target_q   = disco_quantile(mytar, G, q_min, q_max)
    out.target_cdf = cdf_builder(mytar, grid)

    //------------------------------------------------------------------
    // 2. Draw a bootstrap sample for each control in period t
    //------------------------------------------------------------------
    uid  = get_unique(idt)
    cids = select(uid, uid:!=target_id)
    J = length(cids)

    mycon_q   = J(G, J, .)
    mycon_cdf = J(G, J, .)

    for (ci=1; ci<=J; ci++) {
        cd     = select(yt, idt:==cids[ci])
        c_len  = length(cd)
        indices_c = ceil(runiform(c_len,1):*c_len)
        mycon    = cd[indices_c]

        // For this control's bootstrapped sample, compute G-point quantile & CDF
        mycon_q[,ci]   = disco_quantile(mycon, G, q_min, q_max)
        mycon_cdf[,ci] = cdf_builder(mycon, grid)
    }

    //------------------------------------------------------------------
    // 3. If pre-treatment (t <= T0), solve for weights
    //------------------------------------------------------------------
    //    Note: We no longer do "controls_resampled = ..." 
    //    because we now have mycon_q, mycon_cdf already.
    //------------------------------------------------------------------
    out.weights = J(J, 1, .)

    if (t <= T0) {
        if (mixture == 0) {
            // Quantile-based solver with precomputed quantiles
            out.weights = disco_solve_weights(mycon_q, out.target_q, simplex)
        } 
        else {
            // Mixture-based solver with precomputed CDF
            out.weights = disco_mixture_weights(mycon_cdf, out.target_cdf, simplex)
        }
    }

    //------------------------------------------------------------------
    // 4. Store them so higher-level code can build synthetic distribution
    //------------------------------------------------------------------
    out.controls_q   = mycon_q
    out.controls_cdf = mycon_cdf

    return(out)
}


// Compute bootstrap counterfactuals
struct boot_out bootCounterfactuals(struct iter_out vector iter_results,
                                  real scalar T0, T_max, G,
                                  real vector grid,
                                  real scalar mixture) 
{
    struct boot_out scalar bo
    real scalar J, t, gg, m2, gg3, p
    real vector W_avg, target_q, target_cdf, Q_synth, C_synth, Q_synth_sorted
    real matrix weights_all, quantile_diff, cdf_diff, mycon_q, mycon_cdf

    J = cols(iter_results[1].controls_q)

    quantile_diff = J(G,T_max,.)
    cdf_diff = J(G,T_max,.)
    real matrix quantile_synth, cdf_synth, quantile_t, cdf_t
    quantile_synth = J(G,T_max,.)
    cdf_synth = J(G,T_max,.)
    quantile_t = J(G,T_max,.)
    cdf_t = J(G,T_max,.)

    weights_all = J(T0-1,J,.)
    real scalar t_loop
    for (t_loop=1; t_loop<=T0-1; t_loop++) {
        weights_all[t_loop,.] = (iter_results[t_loop].weights)
    }
    W_avg = (colsum(weights_all)/(T0-1))'

    real scalar t_loop2
    for (t_loop2=1; t_loop2<=T_max; t_loop2++) {
        target_q = iter_results[t_loop2].target_q
        target_cdf = iter_results[t_loop2].target_cdf
        mycon_q = iter_results[t_loop2].controls_q
        mycon_cdf = iter_results[t_loop2].controls_cdf

        quantile_t[,t_loop2] = target_q
        cdf_t[,t_loop2] = target_cdf

        if (mixture==0) {
            Q_synth = mycon_q*W_avg
            Q_synth_sorted = sort(Q_synth,1)
            C_synth = J(G,1,0)
            real scalar gg_loop
            for (gg_loop=1; gg_loop<=G; gg_loop++) {
                C_synth[gg_loop] = sum(Q_synth_sorted:<=grid[gg_loop])/G
            }
        } else {
            C_synth = mycon_cdf*W_avg
            Q_synth = J(G,1,.)
            real scalar m2_loop
            for (m2_loop=1; m2_loop<=G; m2_loop++) {
                p = (m2_loop-1)/(G-1)
                gg3 = 1
                while (gg3<=G & C_synth[gg3]<p) gg3++
                Q_synth[m2_loop] = grid[gg3>G ? G : gg3]
            }
        }

        quantile_synth[,t_loop2] = Q_synth
        cdf_synth[,t_loop2] = C_synth
        quantile_diff[,t_loop2] = target_q - Q_synth
        cdf_diff[,t_loop2] = target_cdf - C_synth
    }

    bo.quantile_diff = quantile_diff
    bo.cdf_diff = cdf_diff
    bo.quantile_synth = quantile_synth
    bo.cdf_synth = cdf_synth
    bo.quantile_t = quantile_t
    bo.cdf_t = cdf_t
    return(bo)
}

// Bootstrap CI
struct CI_out scalar disco_bootstrap_CI(real vector y, real vector id, real vector tt,
                                      real scalar target_id, T0, T_max, G,
                                      real scalar q_min, q_max, simplex, mixture, boots, cl, uniform,
									  real matrix quantile_diff, cdf_diff) 
{
    struct CI_out scalar co
    struct boot_out scalar bo
    struct disco_out scalar main_run
    struct iter_out vector iter_results
    real scalar amin, amax, N, b, t, alpha, lower_idx, upper_idx, idx, b2, q_crit, c_crit
    real vector grid, vals, tmp, qmax_abs, cmax_abs
    real matrix quantile_diff_boot, cdf_diff_boot, qdiff_lower, qdiff_upper, cdiff_lower, cdiff_upper
    real matrix quantile_synth_boot, cdf_synth_boot, quantile_t_boot, cdf_t_boot

    iter_results = iter_out(T_max)
    amin = min(y)
    amax = max(y)
    grid = range(amin, amax, (amax - amin)/(G-1))'

    // main_run = disco_full_run(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)

    N = length(y)
    quantile_diff_boot = J(G*T_max, boots, .)
    cdf_diff_boot = J(G*T_max, boots, .)
    quantile_synth_boot = J(G*T_max, boots, .)
    cdf_synth_boot = J(G*T_max, boots, .)
    quantile_t_boot = J(G*T_max, boots, .)
    cdf_t_boot = J(G*T_max, boots, .)

    real scalar b_loop
    for (b_loop=1; b_loop<=boots; b_loop++) {
        real scalar t_loop
        for (t_loop=1; t_loop<=T_max; t_loop++) {
            struct iter_out scalar out
            out = disco_CI_iter(y, id, tt, target_id, t_loop, T0, G, grid, q_min, q_max, simplex, mixture)
            iter_results[t_loop] = out
        }

        bo = bootCounterfactuals(iter_results, T0, T_max, G, grid, mixture)

        quantile_diff_boot[,b_loop] = vec(bo.quantile_diff)
        cdf_diff_boot[,b_loop] = vec(bo.cdf_diff)
        quantile_synth_boot[,b_loop] = vec(bo.quantile_synth)
        cdf_synth_boot[,b_loop] = vec(bo.cdf_synth)
        quantile_t_boot[,b_loop] = vec(bo.quantile_t)
        cdf_t_boot[,b_loop] = vec(bo.cdf_t)
    }

    qdiff_lower = J(G,T_max,.)
    qdiff_upper = J(G,T_max,.)
    cdiff_lower = J(G,T_max,.)
    cdiff_upper = J(G,T_max,.)

    if (uniform==0) {
        real scalar t_loop2, idx_loop3
        for (t_loop2=1; t_loop2<=T_max; t_loop2++) {
            for (idx_loop3=1; idx_loop3<=G; idx_loop3++) {
                vals = quantile_diff_boot[(t_loop2-1)*G + idx_loop3,.]'
                tmp = sort(vals, 1)

                alpha = (1-cl)/2
                lower_idx = max((ceil(alpha*boots), 1))
                upper_idx = min((ceil((1-alpha)*boots), boots))

                qdiff_lower[idx_loop3,t_loop2] = tmp[lower_idx]
                qdiff_upper[idx_loop3,t_loop2] = tmp[upper_idx]
            }

            real scalar idx_loop4
            for (idx_loop4=1; idx_loop4<=G; idx_loop4++) {
                vals = cdf_diff_boot[(t_loop2-1)*G + idx_loop4,.]'
                tmp = sort(vals, 1)

                alpha = (1-cl)/2
                lower_idx = max((ceil(alpha*boots), 1))
                upper_idx = min((ceil((1-alpha)*boots), boots))

                cdiff_lower[idx_loop4,t_loop2] = tmp[lower_idx]
                cdiff_upper[idx_loop4,t_loop2] = tmp[upper_idx]
            }
        }
    } 
    else {

        qmax_abs = J(boots,1,.)
        cmax_abs = J(boots,1,.)

        real scalar b2_loop
        for (b2_loop=1; b2_loop<=boots; b2_loop++) {
            real matrix qdiff_mat, cdiff_mat, qdiff_err, cdiff_err

            qdiff_mat = rowshape(quantile_diff_boot[,b2_loop], G)
            cdiff_mat = rowshape(cdf_diff_boot[,b2_loop], G)

            qdiff_err = qdiff_mat - quantile_diff
            cdiff_err = cdiff_mat - cdf_diff

            qmax_abs[b2_loop] = max(abs(vec(qdiff_err)))
            cmax_abs[b2_loop] = max(abs(vec(cdiff_err)))
        }

        tmp = sort(qmax_abs, 1)
        q_crit = tmp[ceil((1-(1-cl)/2)*boots)]
        tmp = sort(cmax_abs, 1)
        c_crit = tmp[ceil((1-(1-cl)/2)*boots)]

        real scalar t_loop3, m_i2, g_i2, base_val2
        for (t_loop3=1; t_loop3<=T_max; t_loop3++) {
            for (m_i2=1; m_i2<=G; m_i2++) {
                base_val2 = quantile_diff[m_i2,t_loop3]
                qdiff_lower[m_i2,t_loop3] = base_val2 - q_crit
                qdiff_upper[m_i2,t_loop3] = base_val2 + q_crit
            }
            for (g_i2=1; g_i2<=G; g_i2++) {
                base_val2 = cdf_diff[g_i2,t_loop3]
                cdiff_lower[g_i2,t_loop3] = base_val2 - c_crit
                cdiff_upper[g_i2,t_loop3] = base_val2 + c_crit
            }
        }
    }

    co.qdiff_lower = qdiff_lower
    co.qdiff_upper = qdiff_upper
    co.cdiff_lower = cdiff_lower
    co.cdiff_upper = cdiff_upper
	


    return(co)
}

// Wrapper for main DISCO
real scalar disco_wrapper(real vector y, id, tt,
                         real scalar target_id, T0, T_max, G, q_min, q_max, simplex, mixture) 
{
    struct disco_out scalar results
    results = disco_full_run(y, id, tt, target_id, T0, T_max, G, q_min, q_max, simplex, mixture)

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

// Wrapper for CI
real scalar disco_ci_wrapper(real vector y, id, tt,
                           real scalar target_id, T0, T_max, G,
                           real scalar q_min, q_max, simplex, mixture, boots, cl, uniform, 
						   real matrix quantile_diff, real matrix cdf_diff)
{
    struct CI_out scalar results

    results = disco_bootstrap_CI(y, id, tt, target_id, T0, T_max, G,
                               q_min, q_max, simplex, mixture,
                               boots, cl, uniform, quantile_diff, cdf_diff)

    st_matrix("qdiff_lower", results.qdiff_lower)
    st_matrix("qdiff_upper", results.qdiff_upper)
    st_matrix("cdiff_lower", results.cdiff_lower)
    st_matrix("cdiff_upper", results.cdiff_upper)

    return(0)
}



// Compute summary stats
real scalar compute_summary_stats(string scalar agg, real vector sample_points, 
                                real scalar T0, real scalar T_max, real matrix quantile_diff,
                                real matrix cdf_diff, real scalar CI, real scalar cl) 
{
    if (!anyof(("quantile", "cdf", "quantileDiff", "cdfDiff"), agg)) {
        errprintf("Invalid aggregation type\n")
        return(1)
    }

    real scalar is_cdf
    is_cdf = (agg == "cdf" | agg == "cdfDiff")

    real vector grid_points
    if (is_cdf) {
        real scalar amin, amax
        amin = st_numscalar("amin")
        amax = st_numscalar("amax")

        real vector grid_points_temp
        if (max(sample_points) <= 1 & min(sample_points) >= 0) {
            grid_points_temp = amin :+ sample_points :* (amax - amin)
            if (min(grid_points_temp) > amin) grid_points = amin \ grid_points_temp'
            else grid_points = grid_points_temp'
            if (max(grid_points_temp) < amax) grid_points = grid_points \ amax
        } else {
            grid_points = sort(sample_points', 1)
        }
    } else {
        grid_points = sort(sample_points', 1)
        if (min(grid_points) > 0) grid_points = 0 \ grid_points
        if (max(grid_points) < 1) grid_points = grid_points \ 1
    }

    real scalar n_intervals
    n_intervals = length(grid_points) - 1

    real matrix summary_stats
    summary_stats = J(n_intervals * (T_max - T0 + 1), 7, .)

    real scalar row, G, M, t, i
    row = 1

    real vector grid, idx, prob_grid

    for(t = T0; t <= T_max; t++) {
        for(i = 1; i <= n_intervals; i++) {
            summary_stats[row, 1] = t
            summary_stats[row, 2] = grid_points[i]
            summary_stats[row, 3] = grid_points[i + 1]

            if (is_cdf) {
                amin = st_numscalar("amin")
                amax = st_numscalar("amax")
                G = rows(cdf_diff)
                grid = range(amin, amax, (amax - amin)/(G-1))'
                idx = selectindex(grid :>= grid_points[i] :& grid :<= grid_points[i + 1])
            } else {
                M = rows(quantile_diff)
                prob_grid = range(0, 1, 1/(M-1))'
                idx = selectindex(prob_grid :>= grid_points[i] :& prob_grid :<= grid_points[i + 1])
            }

            if (length(idx) > 0) {
                if (is_cdf) {
                    summary_stats[row, 4] = mean(cdf_diff[idx, t])
                } else {
                    summary_stats[row, 4] = mean(quantile_diff[idx, t])
                }
            }

            if (CI) {
                real matrix diff_lower, diff_upper
                if (is_cdf) {
                    diff_lower = st_matrix("cdiff_lower")
                    diff_upper = st_matrix("cdiff_upper")
                } else {
                    diff_lower = st_matrix("qdiff_lower")
                    diff_upper = st_matrix("qdiff_upper")
                }

                if (length(idx) > 0) {
                    if (is_cdf) {
                        summary_stats[row, 4] = mean(cdf_diff[idx, t])
                        summary_stats[row, 6] = mean(diff_lower[idx, t])
                        summary_stats[row, 7] = mean(diff_upper[idx, t])
                        summary_stats[row, 5] = (summary_stats[row, 7] - summary_stats[row, 6])/(2*1.96)
                    } else {
                        summary_stats[row, 4] = mean(quantile_diff[idx, t])
                        summary_stats[row, 6] = mean(diff_lower[idx, t])
                        summary_stats[row, 7] = mean(diff_upper[idx, t])
                        summary_stats[row, 5] = (summary_stats[row, 7] - summary_stats[row, 6])/(2*1.96)
                    }
                }
            }

            row++
        }
    }

    st_matrix("summary_stats", summary_stats)

    return(0)
}
end
