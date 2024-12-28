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
real vector disco_quantile(real vector X, real scalar M, real scalar q_min, real scalar q_max) {
    real vector p
    real scalar j, p_j
    
    p = J(M,1,.)
    real scalar j2
    for (j2=1; j2<=M; j2++) {
        p_j = q_min + (q_max - q_min)*(j2-1)/(M-1)
        p[j2] = p_j
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
    
    p = J(M,1,.)
    real scalar i2
    for (i2=1; i2<=M; i2++) {
        p[i2] = q_min + (q_max - q_min)*(i2-1)/(M-1)
    }

    real matrix Cq
    Cq = J(M,J,.)
    real scalar i3
    for (i3=1; i3<=J; i3++) {
        Cq[,i3] = disco_quantile_points(controls[,i3], p)
    }
    target_s = disco_quantile_points(target, p)

    C_in = Cq
    G = 2*(C_in' * C_in)
    g0 = -2*(C_in' * target_s)

    CE = J(J,1,1)
    ce0 = (-1)

    if (simplex==1) {
        CI = I(J)
        ci0 = J(J,1,0)
    } else {
        CI = J(J,1,1)
        ci0 = (-1e20)
    }

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
    real scalar i, J, n, grid_min, grid_max, ss, val, g, total_vars
    
    real vector alldata, grid_rand, target_cdf, w, sol
    real matrix control_cdf, c, ecmat, lowerbd, upperbd
    real vector bec

    n = length(target)
    J = cols(controls)

    alldata = target
    real scalar i4
    for (i4=1; i4<=J; i4++) {
        alldata = alldata \ controls[,i4]
    }

    grid_min = min(alldata)
    grid_max = max(alldata)
    grid_rand = runiform(M,1)*(grid_max - grid_min) + J(M,1, grid_min)

    target_cdf = cdf_at_points(target, grid_rand)
    control_cdf = J(M,J,.)
    real scalar i5
    for (i5=1; i5<=J; i5++) {
        control_cdf[,i5] = cdf_at_points(controls[,i5], grid_rand)
    }

    total_vars = J + 2*M

    c = J(1, total_vars, 0)
    c[1,(J+1)..(J+M)] = J(1,M,1)
    c[1,(J+M+1)..(J+2*M)] = J(1,M,1)

    ecmat = J(M+1, total_vars, 0)
    bec = J(M+1,1,0)

    ecmat[1,1..J] = J(1,J,1)
    bec[1] = 1

    real scalar g2
    for (g2=1; g2<=M; g2++) {
        ecmat[1+g2,1..J] = control_cdf[g2,.]
        ecmat[1+g2,J+g2] = -1
        ecmat[1+g2,J+M+g2] = 1
        bec[1+g2] = target_cdf[g2]
    }

    lowerbd = J(1,total_vars,.)
    upperbd = J(1,total_vars,.)

    if (simplex==1) {
        lowerbd[1,1..J] = J(1,J,0)
    }
    lowerbd[1,(J+1)..(J+M)] = J(1,M,0)
    lowerbd[1,(J+M+1)..(J+2*M)] = J(1,M,0)

    class LinearProgram scalar q
    q = LinearProgram()
    q.setCoefficients(c)
    q.setMaxOrMin("min")
    q.setEquality(ecmat,bec)
    q.setBounds(lowerbd,upperbd)

    val = q.optimize()
    if (q.errorcode()!=0) {
        w = J(J,1,1/J)
        return(w)
    }

    sol = q.parameters()
    w = sol[1..J]'
    ss = sum(w)
    if (abs(ss-1)>1e-8) w = w/ss

    return(w)
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


// Main DISCO function
struct disco_out disco_full_run(real vector y, real vector id, real vector tt,
                              real scalar target_id, real scalar T0, real scalar T_max, 
                              real scalar M, real scalar G, real scalar q_min, real scalar q_max,
                              real scalar simplex, real scalar mixture) 
{
    real vector uid, cids, yt, idt, target_data, cd, w, w_temp, W_avg, 
             Q_synth, Q_synth_sorted, cids_t, grid, Tq, Tc, C_synth
    
    real matrix controls_data, Cq, Cc, quantile_diff, cdf_diff, weights_store, 
             period_weights, Q_target_all, Q_synth_all, C_target_all, C_synth_all

    real scalar J_sc, amin, amax, gg, pos, t, ci, m, p, ci2, gg2, gg3, m2, J_t, ci5, ci6, ci7, t2

    uid = get_unique(id)
    cids = select(uid, uid:!=target_id)
    J_sc = length(cids)

    amin = min(y)
    amax = max(y)
    grid = range(amin, amax, (amax - amin) / (G - 1))'

    st_numscalar("amin", amin)
    st_numscalar("amax", amax)

    quantile_diff = J(M,T_max,.)
    cdf_diff = J(G,T_max,.)
    weights_store = J(T0-1,J_sc,.)
    period_weights = J(T_max,J_sc,.)
    Q_target_all = J(M,T_max,.)
    Q_synth_all = J(M,T_max,.)
    C_target_all = J(G,T_max,.)
    C_synth_all = J(G,T_max,.)

    real scalar t_loop
    for (t_loop=1; t_loop<=T_max; t_loop++) {
        yt = select(y, tt:==t_loop)
        idt = select(id, tt:==t_loop)
        target_data = select(yt, idt:==target_id)
        controls_data = .

        real scalar ci2_loop
        for (ci2_loop=1; ci2_loop<=J_sc; ci2_loop++) {
            cd = select(yt, idt:==cids[ci2_loop])
            if (ci2_loop==1) {
                controls_data = cd
            }
            else {
                controls_data = controls_data, cd
            }
        }

        Tq = disco_quantile(target_data,M,q_min,q_max)
        Tc = cdf_builder(target_data,grid)

        if (t_loop<=T0-1) {
            if (mixture==0) {
                w = disco_solve_weights(controls_data, target_data, M, q_min, q_max, "", 7, simplex)
            } else {
                w_temp = disco_mixture_weights(controls_data, target_data, M, simplex)
                w = w_temp'
            }
            weights_store[t_loop,.] = w
            period_weights[t_loop,.] = w
        } else {
            period_weights[t_loop,.] = J(J_sc,1,.)'
        }

        Q_target_all[,t_loop] = Tq
        C_target_all[,t_loop] = Tc
    }

    W_avg = (colsum(weights_store)/((T0-1)))'
    real scalar t2_loop
    for (t2_loop=T0; t2_loop<=T_max; t2_loop++) {
        period_weights[t2_loop,.] = W_avg'
    }

    real scalar t_loop2
    for (t_loop2=1; t_loop2<=T_max; t_loop2++) {
        yt = select(y, tt:==t_loop2)
        idt = select(id, tt:==t_loop2)
        target_data = select(yt, idt:==target_id)
        cids_t = select(get_unique(idt), get_unique(idt):!=target_id)

        J_t = length(cids_t)
        controls_data = .
        real scalar ci5_loop
        for (ci5_loop=1; ci5_loop<=J_t; ci5_loop++) {
            cd = select(yt, idt:==cids_t[ci5_loop])
            if (ci5_loop==1) controls_data = cd
            else controls_data = controls_data, cd
        }

        w = period_weights[t_loop2,.]'

        if (mixture==0) {
            real matrix Cq2
            Cq2 = J(M,J_t,.)
            real scalar ci6_loop
            for (ci6_loop=1; ci6_loop<=J_t; ci6_loop++) {
                Cq2[,ci6_loop] = disco_quantile(controls_data[,ci6_loop],M,q_min,q_max)
            }
            Q_synth = Cq2*w
            Q_synth_sorted = sort(Q_synth, 1)

            C_synth = J(G,1,0)
            real scalar gg2_loop
            for (gg2_loop=1; gg2_loop<=G; gg2_loop++) {
                pos = sum(Q_synth_sorted:<=grid[gg2_loop])
                C_synth[gg2_loop] = pos/M
            }
        } else {
            real matrix Cc2
            Cc2 = J(G,J_t,.)
            real scalar ci7_loop
            for (ci7_loop=1; ci7_loop<=J_t; ci7_loop++) {
                Cc2[,ci7_loop] = cdf_builder(controls_data[,ci7_loop],grid)
            }
            C_synth = Cc2*w
            Q_synth = J(M,1,.)

            real scalar m2_loop
            for (m2_loop=1; m2_loop<=M; m2_loop++) {
                p = (m2_loop-1)/(M-1)
                gg3 = 1
                while (gg3<G && C_synth[gg3]<p) {
                    gg3++
                }
                if (gg3>G) gg3=G
                Q_synth[m2_loop] = grid[gg3]
            }
        }

        Q_synth_all[,t_loop2] = Q_synth
        C_synth_all[,t_loop2] = C_synth
        quantile_diff[,t_loop2] = Q_target_all[,t_loop2] - Q_synth
        cdf_diff[,t_loop2] = C_target_all[,t_loop2] - C_synth
    }

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

// Compute ratio
real scalar disco_compute_ratio(real vector y, id, tt, 
                              real scalar target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture) 
{
    struct disco_out scalar rr
    rr = disco_full_run(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)
    
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
                                 real scalar target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture) 
{
    real scalar actual_ratio, pval, rj, J, count, j
    actual_ratio = disco_compute_ratio(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)
    
    real vector uid, cids
    uid = get_unique(id)
    cids = select(uid, uid:!=target_id)
    J = length(cids)

    count = 0
    real scalar j_loop
    for (j_loop=1; j_loop<=J; j_loop++) {
        rj = disco_compute_ratio(y, id, tt, cids[j_loop], T0, T_max, M, G, q_min, q_max, simplex, mixture)
        if (rj>=actual_ratio) count = count + 1
    }

    pval = (count+1)/(J+1)
    return(pval)
}

// Compute CI bounds
struct CI_out scalar compute_CI_bounds(real matrix quantile_diff_boot, cdf_diff_boot,
                                     real scalar M, G, T_max, cl, uniform,
                                     struct disco_out scalar main_run)
{
    struct CI_out scalar co
    co.qdiff_lower = J(M, T_max, .)
    co.qdiff_upper = J(M, T_max, .)
    co.cdiff_lower = J(G, T_max, .)
    co.cdiff_upper = J(G, T_max, .)

    if(!uniform) {
        real scalar alpha, lower_idx, upper_idx, idx
        real vector vals

        alpha = (1-cl)/2
        lower_idx = max((1, ceil(alpha*cols(quantile_diff_boot))))
        upper_idx = min((cols(quantile_diff_boot), ceil((1-alpha)*cols(quantile_diff_boot))))

        real scalar t, idx_loop
        for(t=1; t<=T_max; t++) {
            for(idx_loop=1; idx_loop<=M; idx_loop++) {
                vals = quantile_diff_boot[(t-1)*M + idx_loop,.]'
                vals = sort(vals, 1)
                co.qdiff_lower[idx_loop,t] = vals[lower_idx]
                co.qdiff_upper[idx_loop,t] = vals[upper_idx]
            }

            real scalar idx_loop2
            for(idx_loop2=1; idx_loop2<=G; idx_loop2++) {
                vals = cdf_diff_boot[(t-1)*G + idx_loop2,.]'
                vals = sort(vals, 1)
                co.cdiff_lower[idx_loop2,t] = vals[lower_idx]
                co.cdiff_upper[idx_loop2,t] = vals[upper_idx]
            }
        }
    }
    else {
        real scalar q_crit, c_crit, m_i, g_i, base_val
        // Need loop variables declared here
        real scalar b, t2

        real vector qmax_abs, cmax_abs
        real matrix qdiff_mat, cdiff_mat, qdiff_err, cdiff_err

        qmax_abs = J(cols(quantile_diff_boot), 1, .)
        cmax_abs = J(cols(quantile_diff_boot), 1, .)

        for(b=1; b<=cols(quantile_diff_boot); b++) {
            qdiff_mat = rowshape(quantile_diff_boot[,b], M)
            cdiff_mat = rowshape(cdf_diff_boot[,b], G)
            
            qdiff_err = qdiff_mat :- main_run.quantile_diff
            cdiff_err = cdiff_mat :- main_run.cdf_diff

            qmax_abs[b] = max(abs(vec(qdiff_err)))
            cmax_abs[b] = max(abs(vec(cdiff_err)))
        }

        real vector tmp
        tmp = sort(qmax_abs, 1)
        q_crit = tmp[ceil((1-(1-cl)/2)*cols(quantile_diff_boot))]
        tmp = sort(cmax_abs, 1)
        c_crit = tmp[ceil((1-(1-cl)/2)*cols(quantile_diff_boot))]

        for(t2=1; t2<=T_max; t2++) {
            for(m_i=1; m_i<=M; m_i++) {
                base_val = main_run.quantile_diff[m_i,t2]
                co.qdiff_lower[m_i,t2] = base_val - q_crit
                co.qdiff_upper[m_i,t2] = base_val + q_crit
            }

            for(g_i=1; g_i<=G; g_i++) {
                base_val = main_run.cdf_diff[g_i,t2]
                co.cdiff_lower[g_i,t2] = base_val - c_crit
                co.cdiff_upper[g_i,t2] = base_val + c_crit
            }
        }
    }

    return(co)
}

// Bootstrap iteration
struct iter_out disco_CI_iter(real vector y, real vector id, real vector tt,
                            real scalar target_id, t, T0, M, G,
                            real vector grid,
                            real scalar q_min, q_max, simplex, mixture) 
{
    struct iter_out scalar out
    real scalar t_len, c_len, J
    real scalar ci, cc
    real vector yt, idt, target_data, mytar, uid, cids, cd, mycon, indices_t, indices_c
    real matrix mycon_q, mycon_cdf, controls_resampled

    yt = select(y, tt:==t)
    idt = select(id, tt:==t)
    target_data = select(yt, idt:==target_id)

    t_len = length(target_data)
    indices_t = ceil(runiform(t_len,1):*t_len)
    mytar = target_data[indices_t]

    out.target_q = disco_quantile(mytar, M, q_min, q_max)
    out.target_cdf = cdf_builder(mytar, grid)

    uid = get_unique(idt)
    cids = select(uid, uid:!=target_id)
    J = length(cids)

    mycon_q = J(M, J, .)
    mycon_cdf = J(G, J, .)

    for (ci=1; ci<=J; ci++) {
        cd = select(yt, idt:==cids[ci])
        c_len = length(cd)
        indices_c = ceil(runiform(c_len,1):*c_len)
        mycon = cd[indices_c]

        mycon_q[,ci] = disco_quantile(mycon, M, q_min, q_max)
        mycon_cdf[,ci] = cdf_builder(mycon, grid)
    }

    out.weights = J(J, 1, .)

    if (t<=T0) {
        if (mixture==0) {
            out.weights = disco_solve_weights(mycon, mytar, M, q_min, q_max, "", 7, simplex)

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
            out.weights = disco_solve_weights(controls_resampled, mytar, M, q_min, q_max, "", 7, simplex)
        } 
        else {
            out.weights = disco_mixture_weights(mycon_cdf, out.target_cdf, M, simplex)
        }
    }

    out.controls_q = mycon_q
    out.controls_cdf = mycon_cdf

    return(out)
}

// Compute bootstrap counterfactuals
struct boot_out bootCounterfactuals(struct iter_out vector iter_results,
                                  real scalar T0, T_max, M, G,
                                  real vector grid,
                                  real scalar mixture) 
{
    struct boot_out scalar bo
    real scalar J, t, gg, m2, gg3, p
    real vector W_avg, target_q, target_cdf, Q_synth, C_synth, Q_synth_sorted
    real matrix weights_all, quantile_diff, cdf_diff, mycon_q, mycon_cdf

    J = cols(iter_results[1].controls_q)

    quantile_diff = J(M,T_max,.)
    cdf_diff = J(G,T_max,.)
    real matrix quantile_synth, cdf_synth, quantile_t, cdf_t
    quantile_synth = J(M,T_max,.)
    cdf_synth = J(G,T_max,.)
    quantile_t = J(M,T_max,.)
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
                C_synth[gg_loop] = sum(Q_synth_sorted:<=grid[gg_loop])/M
            }
        } else {
            C_synth = mycon_cdf*W_avg
            Q_synth = J(M,1,.)
            real scalar m2_loop
            for (m2_loop=1; m2_loop<=M; m2_loop++) {
                p = (m2_loop-1)/(M-1)
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
                                      real scalar target_id, T0, T_max, M, G,
                                      real scalar q_min, q_max, simplex, mixture, boots, cl, uniform) 
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

    main_run = disco_full_run(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)

    N = length(y)
    quantile_diff_boot = J(M*T_max, boots, .)
    cdf_diff_boot = J(G*T_max, boots, .)
    quantile_synth_boot = J(M*T_max, boots, .)
    cdf_synth_boot = J(G*T_max, boots, .)
    quantile_t_boot = J(M*T_max, boots, .)
    cdf_t_boot = J(G*T_max, boots, .)

    real scalar b_loop
    for (b_loop=1; b_loop<=boots; b_loop++) {
        real scalar t_loop
        for (t_loop=1; t_loop<=T_max; t_loop++) {
            struct iter_out scalar out
            out = disco_CI_iter(y, id, tt, target_id, t_loop, T0, M, G, grid, q_min, q_max, simplex, mixture)
            iter_results[t_loop] = out
        }

        bo = bootCounterfactuals(iter_results, T0, T_max, M, G, grid, mixture)

        quantile_diff_boot[,b_loop] = vec(bo.quantile_diff)
        cdf_diff_boot[,b_loop] = vec(bo.cdf_diff)
        quantile_synth_boot[,b_loop] = vec(bo.quantile_synth)
        cdf_synth_boot[,b_loop] = vec(bo.cdf_synth)
        quantile_t_boot[,b_loop] = vec(bo.quantile_t)
        cdf_t_boot[,b_loop] = vec(bo.cdf_t)
    }

    qdiff_lower = J(M,T_max,.)
    qdiff_upper = J(M,T_max,.)
    cdiff_lower = J(G,T_max,.)
    cdiff_upper = J(G,T_max,.)

    if (uniform==0) {
        real scalar t_loop2, idx_loop3
        for (t_loop2=1; t_loop2<=T_max; t_loop2++) {
            for (idx_loop3=1; idx_loop3<=M; idx_loop3++) {
                vals = quantile_diff_boot[(t_loop2-1)*M + idx_loop3,.]'
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

            qdiff_mat = rowshape(quantile_diff_boot[,b2_loop], M)
            cdiff_mat = rowshape(cdf_diff_boot[,b2_loop], G)

            qdiff_err = qdiff_mat - main_run.quantile_diff
            cdiff_err = cdiff_mat - main_run.cdf_diff

            qmax_abs[b2_loop] = max(abs(vec(qdiff_err)))
            cmax_abs[b2_loop] = max(abs(vec(cdiff_err)))
        }

        tmp = sort(qmax_abs, 1)
        q_crit = tmp[ceil((1-(1-cl)/2)*boots)]
        tmp = sort(cmax_abs, 1)
        c_crit = tmp[ceil((1-(1-cl)/2)*boots)]

        real scalar t_loop3, m_i2, g_i2, base_val2
        for (t_loop3=1; t_loop3<=T_max; t_loop3++) {
            for (m_i2=1; m_i2<=M; m_i2++) {
                base_val2 = main_run.quantile_diff[m_i2,t_loop3]
                qdiff_lower[m_i2,t_loop3] = base_val2 - q_crit
                qdiff_upper[m_i2,t_loop3] = base_val2 + q_crit
            }
            for (g_i2=1; g_i2<=G; g_i2++) {
                base_val2 = main_run.cdf_diff[g_i2,t_loop3]
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
                         real scalar target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture) 
{
    struct disco_out scalar results
    results = disco_full_run(y, id, tt, target_id, T0, T_max, M, G, q_min, q_max, simplex, mixture)

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

// Wrapper for permutation test
real scalar disco_permutation_wrapper(string scalar y_name, id_name, tt_name,
                                    real scalar target_id, T0, T_max, M, G,
                                    real scalar q_min, q_max, simplex, mixture)
{
    return(0)
}

// Wrapper for CI
real scalar disco_ci_wrapper(real vector y, id, tt,
                           real scalar target_id, T0, T_max, M, G,
                           real scalar q_min, q_max, simplex, mixture, boots, cl, uniform)
{
    struct CI_out scalar results
    results = disco_bootstrap_CI(y, id, tt, target_id, T0, T_max, M, G,
                               q_min, q_max, simplex, mixture,
                               boots, cl, uniform)

    st_matrix("qdiff_lower", results.qdiff_lower)
    st_matrix("qdiff_upper", results.qdiff_upper)
    st_matrix("cdiff_lower", results.cdiff_lower)
    st_matrix("cdiff_upper", results.cdiff_upper)

    return(0)
}

// Helper function
real vector find_grid_indices(real scalar q_start, real scalar q_end, 
                            real scalar grid_length) 
{
    real scalar idx_start, idx_end
    if (q_start < 0 | q_end > 1) {
        _error("Quantile values must be between 0 and 1")
    }

    idx_start = max((1, ceil(q_start * grid_length)))
    idx_end = min((grid_length, floor(q_end * grid_length)))

    return((idx_start \ idx_end))
}

// Compute mean effect
real scalar compute_mean_effect(real matrix data, real scalar idx_start, 
                              real scalar idx_end, real scalar t) 
{
    real vector range
    range = data[|idx_start, t \ idx_end, t|]
    if (anyof(range, .)) {
        return(.)
    }
    return(mean(range))
}

// Compute bootstrap statistics
real vector compute_boot_stats(real matrix boot_data, real scalar idx_start, 
                             real scalar idx_end, real scalar t, 
                             real scalar cl) 
{
    real vector boot_means, stats
    real scalar b
    boot_means = J(cols(boot_data), 1, .)
    for(b=1; b<=cols(boot_data); b++) {
        boot_means[b] = mean(boot_data[|idx_start, t, b \ idx_end, t, b|])
    }

    boot_means = select(boot_means, !missing(boot_means))
    if (length(boot_means) == 0) {
        return(J(3, 1, .))
    }

    real scalar se
    se = sqrt(variance(boot_means))

    real vector ci_bounds
    ci_bounds = get_ci_bounds(boot_means, cl)

    return((se \ ci_bounds))
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

// Compute CI bounds
real vector get_ci_bounds(real vector boot_samples, real scalar cl) 
{
    real vector sorted, bounds
    real scalar n, lower_idx, upper_idx

    sorted = sort(boot_samples, 1)
    n = length(sorted)

    lower_idx = ceil((1-cl)/2 * n)
    upper_idx = ceil((1-(1-cl)/2) * n)

    if (lower_idx < 1) lower_idx = 1
    if (upper_idx > n) upper_idx = n

    bounds = sorted[lower_idx] \ sorted[upper_idx]

    return(bounds)
}

end
