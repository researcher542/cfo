

ds_fast_path(prob, option) = ds_cs_net_get_feasible_path(prob, option, true)

"""
Get a (fast) feasible path with cs_net. 

If fast_flag is true, then return the fastest path in cs_net.

If fast_flag is false, then return a cs_path with minimum charging stops.

We always assume fully charge the battery at the charging station.

return: primal
"""
function ds_cs_net_get_feasible_path(prob::AbsCfoProb, option::DsOption, fast_flag::Bool = true)
    (;src, des, net, ev) = prob
    (;cs_net, cs_vec) = net
    function getdist(cs_net, u, v, lam)::Float64
        edge = get_edge(cs_net, u, v)
        if isnothing(edge)
            return Inf
        else
            if fast_flag
                return edge.min_t
            else
                return 1.0
            end
        end
    end
    cs_path = shortest_path(Astar(), cs_net, src, des; getdist=getdist)
    N = length(cs_path)-2
    # sub_sol_vec = DsSubsolution[]
    sub_sol_vec = [DsSubsolution() for _ in 1:N+1]
    Threads.@threads :greedy for i in 1:N+1
    # @batch for i in 1:N+1
        if cs_path[i] == cs_path[i+1]
            path = [cs_path[i], cs_path[i+1]]
            t_vec = [0.0]
        else
            # path, t_vec, w_cost = ds_recover_primal_edge(prob, cs_path, i, dual, option) 
            src1 = cs_path[i]
            des1 = cs_path[i+1]
            cs_edge = get_edge(cs_net, src1, des1)
            # should use min_t here becase we have refined the edge.
            (path, t_vec, e_cost) = paso(net, src1, des1, cs_edge.min_t; veh=ev, debug_msg = false, output_step = false, visitonce = false, breakearly = false)
            if e_cost > ev.cap
                idx = 2
                while (e_cost > ev.cap) && (idx <= length(cs_edge.t_vec))
                    # @warn "e_cost larger than capacity. Maybe the cs_net is outdated...." src des src1 des1 e_cost ev.cap cs_edge.t_vec[idx-1] cs_edge.c_vec[idx-1]
                    T1 = cs_edge.t_vec[idx]
                    (path, t_vec, e_cost) = paso(net, src1, des1, T1; veh=ev, debug_msg = false, output_step = false, visitonce = false, breakearly = false)
                    idx += 1
                end
            end
        end
        ics_vec = findall(cs->(cs.idx == cs_path[i]), cs_vec)
        @assert(length(ics_vec) == 1)
        ics = ics_vec[1]
        tc = (i==1) ? 0.0 : g_max_charge_time
        tw = (i==1) ? 0.0 : g_min_wait_time
        subsol = DsSubsolution(path, t_vec, tw, tc, 0.0, 0.0, ics)
        sub_sol_vec[i] = subsol
        # push!(sub_sol_vec, subsol)
    end
    primal1 = DsPrimal(;sub_sol_vec=sub_sol_vec)
    primal1 = ds_update_tau_beta!(prob, primal1, option)
    return primal1
end
