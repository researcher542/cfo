
function ds_reopt_cs_net_get_lbub(prob::AbsCfoProb, primal::DsPrimal, option::DsOption)
    (;ev, net) = prob
    (;cs_net) = net
    cs_path = get_cs_path(primal) 

    N = get_N(primal)
    # f = x -> ds_reopt_objective(x, prob, primal, option)
    minmax_t_vec = [minmax_t(cs_net, cs_path[i], cs_path[i+1], ev) for i in 1:N+1]
    t_vec_lb = [t[1] for t in minmax_t_vec]
    t_vec_ub = [t[2] for t in minmax_t_vec]
    tol = 0e-3
    # Note: the length of tc and tw_vec is N+1 because we used vcat(0, ...) here.
    max_wait_time = prob.T
    tc_vec_lb = vcat(0.0, [g_min_charge_time for _ in 1:N])
    tc_vec_ub = vcat(0.0, [g_max_charge_time for _ in 1:N])
    tw_vec_lb = vcat(0.0, [g_min_wait_time for _ in 1:N])
    tw_vec_ub = vcat(0.0, [max_wait_time for _ in 1:N])
    lb = vcat(t_vec_lb, tc_vec_lb, tw_vec_lb)
    ub = vcat(t_vec_ub, tc_vec_ub, tw_vec_ub)
    # @show lb tc_vec_lb tw_vec_lb t_vec_lb
    return lb, ub 
end

"""
Ideally, given any travel time vector, we can directly get tc_vec and tw_vec with minimum objective. 
TODO: this can be tricky, I cannot think of a better way for this...
"""
function ds_reopt_tr_vec2t_vec_all(tr_vec, primal::DsPrimal, prob::AbsCfoProb)
    cs_net = prob.net.cs_net
    B = prob.ev.cap
    beta0 = prob.β0

    cs_path = get_cs_path(primal)
    N = get_N(primal)
    cs_edge_vec = [get_edge(cs_net, cs_path[i], cs_path[i+1]) for i in 1:N+1]
    e_cost_vec =  [edge(tr) for (tr, edge) in zip(tr_vec, cs_edge_vec)]
end


function ds_reopt_cs_net_primal2tvec(primal::DsPrimal)
    N = get_N(primal) 
    sol_vec = primal.sub_sol_vec
    t_vec = [sum(sol.t_vec) for sol in sol_vec]
    tc_vec = [sol.tc for sol in sol_vec]
    tw_vec = [sol.tw for sol in sol_vec]
    return vcat(t_vec, tc_vec, tw_vec) 
end

"""
Decompose the vectorized t vector to sub components.
t_vec_all = [t_vec, tc_vec, tw_vec]
tc_vec, tw_vec have size N+1, where we define tw[1] == tc[1] == 0
"""
function ds_reopt_cs_net_decompose(t_vec_all::Vector{Float64}, primal::DsPrimal)
    N = get_N(primal)
    len_t_vec = N+1
    t_vec = t_vec_all[1:len_t_vec]
    tc_vec = t_vec_all[len_t_vec+1:len_t_vec+N+1]
    tw_vec = t_vec_all[len_t_vec+N+2:end]
    # @show t_vec tc_vec tw_vec
    return t_vec, tc_vec, tw_vec
end

function ds_reopt_cs_net_tvec2primal(t_vec_all, primal::DsPrimal, prob::AbsCfoProb, option)
    (; t_mul, b_mul) = option
    (; net, src, des, ev) = prob
    primal_new = copy(primal)
    t_vec, tc_vec, tw_vec = ds_reopt_cs_net_decompose(t_vec_all, primal) 
    # @debug "recover t,tc,tw from t_vec_all" t_vec tc_vec tw_vec
    # @show t_vec_all
    # @debug "current thread_id is $(Threads.threadid())"
    Threads.@threads :greedy for isol in eachindex(primal_new.sub_sol_vec)
        sol = primal_new.sub_sol_vec[isol]
        ddl::Float64 = t_vec[isol]
        path = sol.path
        if path[1] == path[end]
            path_paso = Int[path[1], path[end]]
            t_vec_paso = Float64[0.0]
        else
            path_paso, t_vec_paso, e_cost = paso(net, sol.path[1], sol.path[end], ddl; veh=ev, output_step = false, debug_msg = false, visitonce = false, breakearly=false)
        end
        sol.path = path_paso
        sol.t_vec = t_vec_paso
        sol.tc = tc_vec[isol]
        sol.tw = tw_vec[isol]
    end

    primal_new = ds_update_tau_beta!(prob, primal_new, option)
    return primal_new 
end

