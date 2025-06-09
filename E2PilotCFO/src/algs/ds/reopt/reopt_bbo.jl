"""
Re-optimize primal with black-box optiomization.
    Note: this function can be slow for large-scale instance
"""
function ds_reoptimize_primal_bbo(prob::AbsCfoProb, primal::DsPrimal, option::DsOption)
    (;net, ev) = prob
    N = get_N(primal)
    f = x -> ds_reopt_objective(x, prob, primal, option)
    minmax_t_vec = reduce(vcat, [
        [ minmax_t(net, path[i], path[i+1], ev) for i in 1:length(path)-1] 
            for path in [sol.path for sol in primal.sub_sol_vec]
        ])
    t_vec_lb = [t[1] for t in minmax_t_vec]
    t_vec_ub = [t[2] for t in minmax_t_vec]
    tc_vec_lb = vcat(0.0, [g_min_charge_time for _ in 1:N])
    tc_vec_ub = vcat(0.0, [g_max_charge_time for _ in 1:N])
    tw_vec_lb = vcat(0.0, [g_min_wait_time for _ in 1:N])
    tw_vec_ub = vcat(0.0, [g_max_wait_time for _ in 1:N])
    lb = vcat(t_vec_lb, tc_vec_lb, tw_vec_lb)
    ub = vcat(t_vec_ub, tc_vec_ub, tw_vec_ub)
    # @show t_vec_lb t_vec_ub

    x0 = (lb+ub)./2

    method = SA[:generating_set_search, :adaptive_de_rand_1_bin_radiuslimited, :adaptive_de_rand_1_bin][2]
    @timeit g_to "reopt_net" res = BBO.bboptimize(f, x0; SearchRange=[(l,u) for (l,u) in zip(lb, ub)], Method = method, TraceMode= :silent, MaxFuncEvals= 100_000) # MaxTime = 1000e-3)

    cost = BBO.best_fitness(res)
    x = BBO.best_candidate(res)

    # f = x->ds_sigma_objective(x, primal, dual, i, ta, prob.start_time, B, option)
    primal_new = ds_reopt_tvec2primal(x, primal, prob, option)
    return cost, primal_new
    
end

"""
Decompose the vectorized t vector to sub components.
t_vec_all = [t_vec, tc_vec, tw_vec]
tc_vec, tw_vec have size N+1, where we define tw[1] == tc[1] == 0
"""
function ds_reopt_decompose(t_vec_all::Vector{Float64}, primal::DsPrimal)
    N = length(primal.sub_sol_vec) - 1
    # path_all = reduce(vcat, [sol.path[1:end-1] for sol in primal.sub_sol_vec])
    # path_all = vcat(path_all, primal.sub_sol_vec[end].path[end])

    # total length of t_vec for travelling.
    len_t_vec = sum([length(sol.t_vec) for sol in primal.sub_sol_vec])
    tc_vec = t_vec_all[len_t_vec+1:len_t_vec+N+1]
    tw_vec = t_vec_all[len_t_vec+N+2:end]
    @assert (length(tc_vec) == N+1) # "$(length(tc_vec)) $N $len_t_vec"

    t_vec_vec = [zeros(length(sol.t_vec)) for sol in primal.sub_sol_vec]
    start_idx = 1
    for (isol,sol) in enumerate(primal.sub_sol_vec)
        n = length(sol.t_vec)
        t_vec_vec[isol] = t_vec_all[start_idx:start_idx+n-1]
        start_idx += n
    end
    @assert  (start_idx-1 == len_t_vec)
    return t_vec_vec, tc_vec, tw_vec
end


function ds_reopt_tvec2primal(t_vec_all, primal::DsPrimal, prob::AbsCfoProb, option)
    (; t_mul, b_mul) = option
    primal_new = copy(primal)
    t_vec, tc_vec, tw_vec = ds_reopt_decompose(t_vec_all, primal) 
    for (isol,sol) in enumerate(primal_new.sub_sol_vec)
        sol.t_vec = t_vec[isol]
        sol.tc = tc_vec[isol]
        sol.tw = tw_vec[isol]
    end
    primal_new = ds_update_tau_beta!(prob, primal_new, option)
    return primal_new 
end

function ds_reopt_objective(t_vec_all, prob::AbsCfoProb, primal::DsPrimal, option::DsOption)
    (;T) = prob
    (;β_lb, t_mul, b_mul) = option
    N = get_N(primal)
    primal_new = ds_reopt_tvec2primal(t_vec_all, primal, prob, option)
    (;obj, beta_vec, tau_vec) = ds_simulate(prob, primal_new; check_flag = false, predict_mode = prob.predict_mode)
    penalty = 1e10
    obj_dual = obj

    # neg(v) = [x < 0 ?  x : 0.0 for x in v]
    # @show sol_beta_vec = [sol.β for sol in primal_new.sub_sol_vec] 
    beta_err_vec = [pos(β_lb - sol.β) for sol in primal_new.sub_sol_vec] * b_mul
    obj_dual += penalty * sum(beta_err_vec)
    obj_dual += penalty * pos(β_lb*b_mul - beta_vec[end]) 
    obj_dual += penalty * sum([x < 0.0 ? -x : 0.0 for x in beta_vec])
    obj_dual += penalty * pos(tau_vec[end] - T*t_mul)
    return obj_dual 
end

