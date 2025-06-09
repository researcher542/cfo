

"""
An algorithm to mimic the behavior of the human practice.

return: (primal, infeasible_flag)
"""
function cfo_practice(prob::AbsCfoProb, tol_ratio::Float64=0.25; verbose::Int = 0)
    (;net, ev, src, des, β0, T) = prob
    cs_vec = net.cs_vec

    function find_shortest_path0(n1, n2)
        function getdist(net, u, v, lambda)
            dis = ep.distance(net, u, v)
            min_t, max_t = minmax_t(net, u, v, ev)
            return min_t
        end
        path = shortest_path(Astar(), net, n1, n2; getdist=getdist, breakearly=false, visitonce=false)
        time_v = [get_minmax_t(net, path[i], path[i+1], ev)[1] for i in 1:length(path)-1]
        return path, time_v
    end

    
    sub_sol_vec = DsSubsolution[]
    battery_tol = tol_ratio * ev.cap

    cur_node = src
    cur_beta = β0
    istage = 1
    used_cs_vec = Int[]
    while istage < 110
        if cur_node == des
            break
        end
        if verbose > 0
            @debug "" istage cur_node des
        end
        path0, t_vec0 = find_shortest_path0(cur_node, des)
        path = Int[cur_node]
        t_vec = Float64[]
        for ipath in 2:length(path0)
            next_node = path0[ipath]
            # way = getway(net, cur_node, next_node)
            t = t_vec0[ipath-1]
            e = energy_cost_on_road(net, cur_node, next_node, ev, t)
            # @debug "" next_node e

            ### detour to charge
            if ((cur_beta - e) < battery_tol) && (next_node != des)
                if verbose > 0
                    @debug istage "detour to charge" cur_node next_node e cur_beta
                end
                cs = find_station_to_charge(net, cur_node, des, used_cs_vec)
                if cs.idx == cur_node
                    path1 = []
                    t_vec1 = []
                else
                    path1, t_vec1 = find_shortest_path0(cur_node, cs.idx)
                end
                path = vcat(path, path1[2:end])
                t_vec = vcat(t_vec, t_vec1)
                
                tc = (istage==1) ? 0.0 : g_max_charge_time
                tw = (istage==1) ? 0.0 : g_min_wait_time
                # ics_vec = findall(cs->(cs.idx == cs_path[i]), cs_vec)
                # @assert(length(ics_vec) == 1)
                # ics = ics_vec[1]
                subsol = DsSubsolution(path, t_vec, tw, tc, 0.0, 0.0, -1)
                push!(sub_sol_vec, subsol)
                istage += 1
                if verbose > 0
                    @debug "route from $(cur_node) charge at cs=$(cs.idx) ipath=$(ipath) istage=$(istage)"
                    @debug "previous path is $path" 
                end

                push!(used_cs_vec, cs.idx)
                cur_beta = ev.cap
                cur_node = cs.idx
                break ## stop iterating the path0
            end # end of if (cur_beta - e) < battery_tol

            cur_beta = cur_beta - e
            push!(path, next_node)
            push!(t_vec, t)
            if next_node == des
                tc = (istage==1) ? 0.0 : g_max_charge_time
                tw = (istage==1) ? 0.0 : g_min_wait_time
                subsol = DsSubsolution(path, t_vec, tw, tc, 0.0, 0.0, -1)
                push!(sub_sol_vec, subsol)
                istage += 1
            end
            cur_node = next_node

        end # for ipath in 2:length(path0)

        if verbose > 0
            @debug "" cur_node istage (path[1], path[end]) cur_beta ev.cap battery_tol
        end
    end
    if istage > 100
        @warn "There seems circular route. Need to check it later." src des T
        @assert false
    end

    primal1 = DsPrimal(;sub_sol_vec=sub_sol_vec)
    primal1 = ds_update_tau_beta!(prob, primal1, DsOption())

    res0 = ds_simulate(prob, primal1; check_flag = false, predict_mode = prob.predict_mode)

    beta_vec = res0.beta_vec
    min_beta = minimum(beta_vec) * 3.6e6 / ev.cap
    if (min_beta < 0)
        tol_ratio += 0.05
        @debug "The result is not battery feasible $(min_beta * 100) %. Increase the tolerance and try again." tol_ratio 
        if (tol_ratio >= 0.8)
            @warn "tol_ratio is too large. The result is not battery feasible. Return the result without checking the battery feasibility."
            return primal1, true
        end
        return cfo_practice(prob, tol_ratio; verbose=verbose)
    end
    
    return primal1, false

end

"""
A heuristic to find the charging station to charge
"""
function find_station_to_charge(net::Network, src::Int, des::Int, used_cs_vec::Vector{Int})
    # src_node = getnode(net, src)
    # des_node = getnode(net, des)
    # max_dis = j2kwh(β0) / kwh_per_km * 1e3
    # cs_vec = deepcopy(net.cs_vec)
    cs_vec = net.cs_vec
    function sortby(cs)
        if cs.idx in used_cs_vec
            return Inf
        end

        dis1 = distance3d(net, cs.idx, src)

        nd_s = getnode(net, src)
        nd_d = getnode(net, des)
        nd_m = getnode(net, cs.idx)
        ##
        angle1 = angle(nd_s, nd_d, nd_m)
        # @show angle1
        ## add panelty if the charging station is not in the direction towards the destination.
        dis2 = angle1 > pi * 0.5 ? 0.0 : 500e3

        return dis1 + dis2
    end
    cs_vec = sort(cs_vec, by=sortby)
    return cs_vec[1]

    # for cs in net.cs_vec
    #     node = getnode(net, cs.idx)
    #     dis = distance3d(node, src_node)
    #     dis2des = distance3d(node, des_node)
    #     # The angle, close to zero is better.
    #     # ang = angle(src_node, des_node, node) / pi * 180
    #     # degree of 30
    #     # ang_flag = ang <= 30
    #     score = -dis + dis2des
    # end

    # return cs_can_vec[idx][1]
end