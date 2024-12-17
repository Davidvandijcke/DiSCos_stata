mata:

real vector unique(real vector x) {
    real vector y
    y = uniqrows(sort(x,1))
    return(y)
}
real vector disco_quantile(real vector X, real scalar M, real scalar q_min, real scalar q_max) {
    // Implements R's type=7 quantile
    real vector Xs, out
    real scalar N, j, p, alpha, floor_alpha, gamma, idx
    
    // Return vector of missing values if X is empty
    if (length(X)==0) {
        return(J(M,1,.))
    }
    
    Xs = sort(X, 1)
    N = length(Xs)
    out = J(M,1,.)

    for (j=1;j<=M;j++) {
        p = q_min + (q_max - q_min)*(j-1)/(M-1)
        if (p<=0) {
            out[j] = Xs[1]
        }
        else if (p>=1) {
            out[j] = Xs[N]
        }
        else {
            alpha = (N-1)*p+1
            floor_alpha = floor(alpha)
            gamma = alpha - floor_alpha
            idx = floor_alpha
            if (idx<1) idx=1
            if (idx>=N) {
                out[j]=Xs[N]
            }
            else {
                out[j] = Xs[idx]*(1-gamma) + Xs[idx+1]*gamma
            }
        }
    }
    return(out)
}

real vector disco_solve_weights(real matrix C, real vector T, real scalar simplex) {
    real vector W
    real scalar s
    W = invsym(C' * C)*C'*T
    if (simplex==1) {
        W = max(W,0)
        s = sum(W)
        if (s==0) {
            W = J(cols(C),1,1/cols(C))
        } else {
            W = W/s
        }
    } else {
        s = sum(W)
        if (abs(s-1)>1e-12) {
            W = W/s
        }
    }
    return(W)
}

real scalar disco_mixture_eval(real vector b) {
    return(0)
}

real vector disco_mixture_weights(real matrix control_cdf, real vector target_cdf, real scalar simplex) {
    real scalar G, J, total_vars, res, ss
    real vector w, c
    real matrix A, A1
    real vector sol

    G = rows(control_cdf)
    J = cols(control_cdf)
    total_vars = J + 2*G

    o = optimize_init()
    optimize_init_evaluator(o, &disco_mixture_eval())
    optimize_init_type(o,"lp")
    optimize_init_nvar(o,total_vars)

    A1 = J(total_vars,1,0)
    for (j=1;j<=J;j++) A1[j]=1
    optimize_add_equality(o,A1,1)

    if (simplex==1) {
        for (j=1;j<=J;j++) {
            optimize_add_lb(o,j,0)
        }
    }

    for (v=(J+1); v<=total_vars; v++) {
        optimize_add_lb(o,v,0)
    }

    for (g=1; g<=G; g++) {
        A = J(total_vars,1,0)
        for (j=1;j<=J;j++) {
            A[j] = control_cdf[g,j]
        }
        A[J+g] = -1
        A[J+G+g]=1
        optimize_add_equality(o,A,target_cdf[g])
    }

    c = J(total_vars,1,0)
    for (g=1; g<=G; g++) {
        c[J+g]=1
        c[J+G+g]=1
    }
    optimize_init_objective(o,c)
    res = optimize(o)
    if (res!=0) {
        return(J(J,1,1/J))
    }
    sol = optimize_solution(o)
    w = sol[1..J]
    ss = sum(w)
    if (abs(ss-1)>1e-8) w=w/ss
    return(w)
}

real vector cdf_builder(real vector X, real vector grid) {
    real vector Xs, out
    real scalar N, G, g, pos
    Xs = sort(X, 1)
    N = length(Xs)
    G = length(grid)
    out = J(G,1,0)
    for (g=1; g<=G; g++) {
        pos = sum(select(Xs, Xs:<=grid[g]))
        out[g]=pos/N
    }
    return(out)
}

struct disco_out {
    real matrix weights
    real matrix quantile_diff
    real matrix cdf_diff
}

struct disco_out disco_full_run(real vector y, real vector id, real vector tt, real scalar target_id, real scalar T0, real scalar T_max, real scalar M, real scalar G, real scalar q_min, real scalar q_max, real scalar simplex, real scalar mixture) {

    real vector uid, cids, yt, idt, target_data, cd, w, W_avg, Q_synth, Q_synth_sorted, cids_t
    real matrix controls_data, Cq, Cc, quantile_diff, cdf_diff, weights_store, period_weights, Q_target_all, Q_synth_all, C_target_all, C_synth_all
    real scalar J, amin, amax, gg, pos, t, ci, m, p
    real vector grid, Tq, Tc, C_synth

    uid = unique(id)
    cids = select(uid, uid:!=target_id)
    J = length(cids)
	


    amin = min(y)
    amax = max(y)
	grid = range(amin, amax, (amax - amin) / (G - 1))'

    quantile_diff = J(M,T_max,.)
    cdf_diff = J(G,T_max,.)
    weights_store = J(T0-1,J,.)
    period_weights = J(T_max,J,.)

    Q_target_all = J(M,T_max,.)
    Q_synth_all = J(M,T_max,.)
    C_target_all = J(G,T_max,.)
    C_synth_all = J(G,T_max,.)

    for (t=1; t<=T_max; t++) {
        yt = select(y, tt:==t)
        idt = select(id, tt:==t)

        
        target_data = select(yt, idt:==target_id)
        controls_data = .
        
        for (ci=1; ci<=J; ci++) {
            cd = select(yt, idt:==cids[ci])
            
            if (ci==1) {
                controls_data = cd
            }
            else {
                controls_data = controls_data, cd
            }
        }

        Tq = disco_quantile(target_data,M,q_min,q_max)
        Tc = cdf_builder(target_data,grid)



        if (t<=T0-1) {
            if (mixture==0) {

                Cq = J(M,J,.)

                for (ci=1; ci<=J; ci++) {
                    Cq[,ci] = disco_quantile(controls_data[,ci],M,q_min,q_max)
                }
                
                w = disco_solve_weights(Cq,Tq,simplex)
            } else {
                Cc = J(G,J,.)
                for (ci=1; ci<=J; ci++) {
                    Cc[,ci] = cdf_builder(controls_data[,ci],grid)
                }
                w = disco_mixture_weights(Cc,Tc,simplex)
            }
            
            weights_store[t,.] = w'
            period_weights[t,.] = w'
        } else {

            period_weights[t,.] = J(rows(w'), cols(w'),.)
        }

        Q_target_all[,t] = Tq
        C_target_all[,t] = Tc
    }

    W_avg = (colsum(weights_store)/((T0-1)))'
    for (t=T0; t<=T_max; t++) {
        period_weights[t,.] = W_avg'
    }

    for (t=1; t<=T_max; t++) {
        yt = select(y, tt:==t)
        idt = select(id, tt:==t)
        target_data = select(yt, idt:==target_id)
        cids_t = select(unique(idt), unique(idt):!=target_id)
        J_t = length(cids_t)
        controls_data = .
        for (ci=1; ci<=J_t; ci++) {
            cd = select(yt, idt:==cids_t[ci])
            if (ci==1) controls_data = cd
            else controls_data = controls_data, cd
        }

        w = period_weights[t,.]'
		
        if (mixture==0) {
            Cq = J(M,J,.)
            for (ci=1; ci<=J_t; ci++) {
                Cq[,ci] = disco_quantile(controls_data[,ci],M,q_min,q_max)
            }
            Q_synth = Cq*w
            Q_synth_sorted = sort(Q_synth, 1)
            C_synth = J(G,1,0)
            for (gg=1; gg<=G; gg++) {
                pos = sum(select(Q_synth_sorted, Q_synth_sorted:<=grid[gg]))
                C_synth[gg] = pos/M
            }
        } else {
            Cc = J(G,J_t,.)
            for (ci=1; ci<=J_t; ci++) {
                Cc[,ci] = cdf_builder(controls_data[,ci],grid)
            }
            C_synth = Cc*w
            Q_synth = J(M,1,.)
            for (m=1;m<=M;m++) {
                p = (m-1)/(M-1)
                gg = 1
                while (gg<=G & C_synth[gg]<p) gg++
                if (gg>G) gg = G
                Q_synth[m] = grid[gg]
            }
        }


        Q_synth_all[,t] = Q_synth
        C_synth_all[,t] = C_synth

        quantile_diff[,t] = Q_target_all[,t] - Q_synth_all[,t]
        cdf_diff[,t] = C_target_all[,t] - C_synth_all[,t]
    }
	// printf("Dimensions of W_avg: %g x %g\n", rows(W_avg), cols(W_avg))


    struct disco_out scalar r
    r.weights = W_avg

    r.quantile_diff = quantile_diff
    r.cdf_diff = cdf_diff
    return(r)
}

real scalar disco_compute_ratio(real vector y, real vector id, real vector tt, real scalar target_id, real scalar T0, real scalar T_max, real scalar M, real scalar G, real scalar q_min, real scalar q_max, real scalar simplex, real scalar mixture) {
    real scalar dist_t, pre_dist, pre_count, post_dist, post_count, ratio
    real matrix quantile_diff
    struct disco_out rr
    rr = disco_full_run(y,id,tt,target_id,T0,T_max,M,G,q_min,q_max,simplex,mixture)
    pre_dist = 0
    pre_count = 0
    post_dist = 0
    post_count = 0
    for (t=1;t<=T_max;t++) {
        dist_t = mean((rr.quantile_diff[,t]:^2))
        if (t<T0) {
            pre_dist = pre_dist + dist_t
            pre_count = pre_count + 1
        }
        if (t>=T0) {
            post_dist = post_dist + dist_t
            post_count = post_count + 1
        }
    }
    if (pre_count==0 | post_count==0) return(.)
    ratio = sqrt((post_dist/post_count))/sqrt((pre_dist/pre_count))
    return(ratio)
}

real scalar disco_permutation_test(real vector y, real vector id, real vector tt, real scalar target_id, real scalar T0, real scalar T_max, real scalar M, real scalar G, real scalar q_min, real scalar q_max, real scalar simplex, real scalar mixture) {
    real scalar actual_ratio, pval, rj
    real vector uid, cids
    real scalar J, count, j
    actual_ratio = disco_compute_ratio(y,id,tt,target_id,T0,T_max,M,G,q_min,q_max,simplex,mixture)
    uid = unique(id)
    cids = select(uid, uid:!=target_id)
    J = length(cids)
    count = 0
    for (j=1;j<=J;j++) {
        rj = disco_compute_ratio(y,id,tt,cids[j],T0,T_max,M,G,q_min,q_max,simplex,mixture)
        if (rj>=actual_ratio) {
			count = count + 1
		}
    }
    pval = (count+1)/(J+1)
    return(pval)
}
 
 

struct CI_out {
	real matrix qdiff_lower
	real matrix qdiff_upper
	real matrix cdiff_lower
	real matrix cdiff_upper
}

// New helper struct for CI iteration
struct iter_out {
	real vector target_q
	real vector target_cdf
	real matrix controls_q
	real matrix controls_cdf
	real vector weights
}

struct boot_out {
	real matrix quantile_diff
	real matrix cdf_diff
}

// Mimic DiSCo_CI_iter in R
struct iter_out disco_CI_iter(real vector y, real vector id, real vector tt, real scalar target_id, real scalar t,
							  real scalar T0, real scalar M, real scalar G, real vector grid,
							  real scalar q_min, real scalar q_max, real scalar simplex, real scalar mixture) {
	struct iter_out scalar out

	yt = select(y, tt:==t)
	idt = select(id, tt:==t)

	target_data = select(yt, idt:==target_id)
	t_len = length(target_data)
	indices_t = ceil(runiform(t_len,1)*t_len)
	mytar = target_data[indices_t]

	out.target_q = disco_quantile(mytar,M,q_min,q_max)
	out.target_cdf = cdf_builder(mytar,grid)

	uid = unique(idt)
	cids = select(uid, uid:!=target_id)
	J = length(cids)

	mycon_q = J(M,J,.)
	mycon_cdf = J(G,J,.)

	for (ci=1; ci<=J; ci++) {
		cd = select(yt, idt:==cids[ci])
		c_len = length(cd)
		indices_c = ceil(runiform(c_len,1)*c_len)
		mycon = cd[indices_c]
		mycon_q[,ci] = disco_quantile(mycon,M,q_min,q_max)
		mycon_cdf[,ci] = cdf_builder(mycon,grid)
	}

	if (t<=T0) {
		if (mixture==0) {
			out.weights = disco_solve_weights(mycon_q,out.target_q,simplex)
		} else {
			out.weights = disco_mixture_weights(mycon_cdf,out.target_cdf,simplex)
		}
	} else {
		out.weights = J(J,1,.)
	}

	out.controls_q = mycon_q
	out.controls_cdf = mycon_cdf
	return(out)
}

// Mimic bootCounterfactuals in R
struct boot_out bootCounterfactuals(struct iter_out vector iter_results, real scalar T0, real scalar T_max,
									real scalar M, real scalar G, real vector grid, real scalar mixture) {
	struct boot_out scalar bo

	// Compute average weights over pre-treatment
	// iter_results is T_max x 1
	J = cols(iter_results[1,1].controls_q)
	weights_all = J(T0-1,J,.)
	for (t=1; t<=T0-1; t++) {
		weights_all[t,.] = (iter_results[t].weights)'
	}
	W_avg = (colsum(weights_all)/(T0-1))'

	quantile_diff = J(M,T_max,.)
	cdf_diff = J(G,T_max,.)

	for (t=1; t<=T_max; t++) {
		target_q = iter_results[t].target_q
		target_cdf = iter_results[t].target_cdf
		mycon_q = iter_results[t].controls_q
		mycon_cdf = iter_results[t].controls_cdf

		if (mixture==0) {
			Q_synth = mycon_q*W_avg
			Q_synth_sorted = sort(Q_synth,1)
			C_synth = J(G,1,0)
			for (gg=1; gg<=G; gg++) {
				pos = sum(select(Q_synth_sorted, Q_synth_sorted:<=grid[gg]))
				C_synth[gg] = pos/M
			}
		} else {
			C_synth = mycon_cdf*W_avg
			Q_synth = J(M,1,.)
			for (m=1;m<=M;m++) {
				p = (m-1)/(M-1)
				gg=1
				while (gg<=G & C_synth[gg]<p) gg++
				if (gg>G) gg=G
				Q_synth[m]=grid[gg]
			}
		}

		quantile_diff[,t] = target_q - Q_synth
		cdf_diff[,t] = target_cdf - C_synth
	}

	bo.quantile_diff = quantile_diff
	bo.cdf_diff = cdf_diff
	return(bo)
}

struct CI_out scalar disco_bootstrap_CI(real vector y, real vector id, real vector tt,
										  real scalar target_id, real scalar T0, real scalar T_max,
										  real scalar M, real scalar G, real scalar q_min, real scalar q_max,
										  real scalar simplex, real scalar mixture, real scalar boots,
										  real scalar cl, real scalar uniform) {
	struct CI_out scalar co
	struct boot_out scalar bo
	struct disco_out scalar main_run  
	struct iter_out vector iter_results
	iter_results = iter_out(T_max)

	amin = min(y)
	amax = max(y)
	grid = range(amin, amax, (amax - amin)/(G-1))'

	// main run without bootstrap
	main_run = disco_full_run(y,id,tt,target_id,T0,T_max,M,G,q_min,q_max,simplex,mixture)

	quantile_diff_boot = J(M*T_max,boots,.)
	cdf_diff_boot = J(G*T_max,boots,.)
	N = length(y)

	for (b=1;b<=boots;b++) {
		// For each period, do a CI iteration
	
		for (t=1; t<=T_max; t++) {
			struct iter_out scalar out

			// Resample at period t
			out = disco_CI_iter(y,id,tt,target_id,t,T0,M,G,grid,q_min,q_max,simplex,mixture)

			iter_results[t] = out

		}

		// Build bootstrapped counterfactuals
		bo = bootCounterfactuals(iter_results,T0,T_max,M,G,grid,mixture)
		quantile_diff_boot[,b] = vec(bo.quantile_diff)
		cdf_diff_boot[,b] = vec(bo.cdf_diff)
	}
	


	qdiff_lower = J(M,T_max,.)
	qdiff_upper = J(M,T_max,.)
	cdiff_lower = J(G,T_max,.)
	cdiff_upper = J(G,T_max,.)

	if (uniform==0) {
		for (t=1;t<=T_max;t++) {
			for (m_i=1;m_i<=M;m_i++) {
				idx = (t-1)*M + m_i
				vals = quantile_diff_boot[idx,.]'
				tmp = sort(vals, 1)
				alpha=(1-cl)/2
				lower_idx=ceil(alpha*boots)
				upper_idx=ceil((1-alpha)*boots)
				if (lower_idx<1) lower_idx=1
				if (upper_idx>boots) upper_idx=boots
				qdiff_lower[m_i,t]=tmp[lower_idx]
				qdiff_upper[m_i,t]=tmp[upper_idx]
			}
			for (g_i=1;g_i<=G;g_i++) {
				idx = (t-1)*G + g_i
				vals = cdf_diff_boot[idx,.]'
				tmp = sort(vals, 1)
				alpha=(1-cl)/2
				lower_idx=ceil(alpha*boots)
				upper_idx=ceil((1-alpha)*boots)
				if (lower_idx<1) lower_idx=1
				if (upper_idx>boots) upper_idx=boots
				cdiff_lower[g_i,t]=tmp[lower_idx]
				cdiff_upper[g_i,t]=tmp[upper_idx]
			}
		}
	} else {
		// Uniform CIs
		qmax_abs = J(boots,1,.)
		cmax_abs = J(boots,1,.)
		for (b=1;b<=boots;b++) {
			qdiff_err = (rowshape(quantile_diff_boot[,b], M) :- main_run.quantile_diff)
			cdiff_err = (rowshape(cdf_diff_boot[,b], G) :- main_run.cdf_diff)
			qmax_abs[b] = max(abs(vec(qdiff_err)))
			cmax_abs[b] = max(abs(vec(cdiff_err)))
		}

		tmp = sort(qmax_abs, 1)
		q_crit = tmp[ceil((1-(1-cl)/2)*boots)]
		tmp = sort(cmax_abs, 1)
		c_crit = tmp[ceil((1-(1-cl)/2)*boots)]

		for (t=1;t<=T_max;t++) {
			for (m_i=1;m_i<=M;m_i++) {
				base_val = main_run.quantile_diff[m_i,t]
				qdiff_lower[m_i,t]=base_val - q_crit
				qdiff_upper[m_i,t]=base_val + q_crit
			}
			for (g_i=1;g_i<=G;g_i++) {
				base_val = main_run.cdf_diff[g_i,t]
				cdiff_lower[g_i,t]=base_val - c_crit
				cdiff_upper[g_i,t]=base_val + c_crit
			}
		}
	}

	co.qdiff_lower = qdiff_lower
	co.qdiff_upper = qdiff_upper
	co.cdiff_lower = cdiff_lower
	co.cdiff_upper = cdiff_upper
	
	return(co)
}

end
