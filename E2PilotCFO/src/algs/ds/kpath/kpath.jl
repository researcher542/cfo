
@with_kw struct KpOption
    K::Int = 16
    restrict_N::Bool = true
    fast_flag::Bool = true
    thread_flag::Bool = true
end

function kpath_get_primal_vec(prob, option, kp_option)
    (;net, src, des, cs_dict, ev) = prob
    (;cs_net) = net
    (;K, restrict_N, fast_flag, thread_flag) = kp_option
    cs_net_src = cs_net.id2idx[src]
    cs_net_des = cs_net.id2idx[des]
    cs_nv = length(cs_net.data)
    distmx = fill(Inf, cs_nv, cs_nv)
    function get_cs_edge_cost(cs_edge, type::ObjCarbon)
        cs_idx = cs_edge.des
        cs = cs_dict[cs_idx]
        c_cost = 2 * cs.avg_price 

        e_cost =  cs_edge.min_cost / ev.cap
        remain_e = ev.cap - cs_edge.min_cost - ev.cap * 0.001
        if remain_e < option.β_lb
            e_cost += 100.0
        end

        ## add penalty for taking more charging stops.
        cost = 0.5 + c_cost + e_cost
        # cost = cs_edge.min_t * cs.avg_price + 3600.0
        return cost
    end
    function get_cs_edge_cost(cs_edge, type::ObjEnergy)
        e_cost =  cs_edge.min_cost / ev.cap
        remain_e = ev.cap - cs_edge.min_cost - ev.cap * 0.001
        if remain_e < option.β_lb
            e_cost += 100.0
        end
        cost = 0.5 + e_cost
        # return get_cs_edge_cost(cs_edge, ObjCarbon())
        return cost
    end

    function get_cs_edge_cost(cs_edge, type::ObjTime)
        cost = cs_edge.min_t
        cost += 3600.0 * 5.0
        return cost
    end


    for (u, edge_vec1) in cs_net.data
        for edge_p in edge_vec1
            v = edge_p[1]
            cs_edge::CsEdge{Float64} = edge_p[2]
            if fast_flag
                distmx[cs_net.id2idx[u], cs_net.id2idx[v]] = cs_edge.min_t
            else
                distmx[cs_net.id2idx[u], cs_net.id2idx[v]] = get_cs_edge_cost(cs_edge, prob.objtype)
            end
            ## try to minimize the number of stops
        end
    end

    yen_state = Graphs.yen_k_shortest_paths(cs_net.g, cs_net_src, cs_net_des, distmx, K) 
    path_vec = yen_state.paths
    primal_vec = DsPrimal[]
    for path in path_vec
        sub_sol_vec = [DsSubsolution(;path=[cs_net.idx2id[path[ind]], cs_net.idx2id[path[ind+1]]]) for ind in 1:length(path)-1]
        # push!(sub_sol_vec, DsSubsolution(;path=))
        primal = DsPrimal(;sub_sol_vec=sub_sol_vec)
        # @show length(get_cs_path(primal)), option.N
        if (option.N != -1) && length(get_cs_path(primal)) > option.N + 1
            if restrict_N
                continue
            elseif length(get_cs_path(primal)) > 2*(option.N+1)
                # Even if we do not restrict N, we still do not want it to be too large.
                continue
            end
        end
        push!(primal_vec, primal)
    end

    @debug "Getting k-shortest path with after_len=$(length(primal_vec))"

    return primal_vec
end

"""
Use k-shortest path on the cs_net to enumerate possible cs_path, then reoptimize them.
restrict_N: If we want to restrcit the number of charging stops.
fast_flag: If we want to simply use the k-fastest path. If false, we will use some herustics to find better path candidates.
"""
function cfo_kpath(prob::AbsCfoProb, option0::DsOption, kp_option::KpOption)
    @debug "calling cfo_kpath"
    option = copy(option0)

    ds_state = DsState()
    ds_preprocess!(prob, ds_state, option)
    
    primal_vec = kpath_get_primal_vec(prob, option, kp_option)
    @debug "reoptimizing with num of paths=$(length(primal_vec))"
    func = (pri) -> ds_reopt(prob, pri, option)
    if kp_option.thread_flag
        sol_candidates = e2map(func, primal_vec)
    else
        sol_candidates = map(func, primal_vec)
    end
    obj, idx = findmin(x->x[1], sol_candidates)
    obj, best_primal = sol_candidates[idx]

    ds_state.obj = obj
    ds_state.obj_f = obj
    ds_state.primal_f = best_primal

    return [ds_state]
end