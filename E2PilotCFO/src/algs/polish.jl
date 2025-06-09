
"""
Find the closest charging station in the direction of destination
"""
function closest_cs(net::Network, i::Int, des::Int, max_dis::Real)
    min_dis = Inf
    idx = -1
    # node = getnode(net, i)
    # We should not go back to charge...
    dis3 = distance3d(net, i, des)
    for cs in net.cs_vec
        dis1 = distance3d(net, i, cs.idx)
        # if dis1 > max_dis && (idx != -1)
        #     continue
        # end
        dis2 = distance3d(net, des, cs.idx)
        # the charging station should be in the middle of i and des
        # if (dis1^2 + dis2^2 > dis3^2) && (idx != -1)
        #     continue
        # end
        dis = dis1 # + dis2/1e3
        #distance3d(net, )
        if dis < min_dis
            idx = cs.idx
            min_dis = dis 
        end
    end
    return idx
end

"""
Find the closest_node in a path
"""
function ep.closest_node(net::Network, i::Int, path::Vector{Int})
    dis_vec = (distance3d(net, i, inode) for inode in path)
    val,idx = findmin(dis_vec) 
    return idx
end

"""
Given path/speed planning, polish the solution to feasible SoC 
"""
function polish(net::Network, path::Vector, t_vec::Vector, ev::EV, β0, start_time::DateTime; T = Inf, γ::Real = 0.2, slow_flag=true)

    n_node = length(path)
    β_vec::Vector{Float64} = [β0]
    τ_vec::Vector{Float64} = [0.0]
    t_vec_new::Vector{Float64} = []
    path_new::Vector{Int} = []
    obj::Float64 = 0.0
    B = ev.cap
    inode = 1 
    des = path[end]
    remain_kwh = j2kwh(γ*B)
    # on average, 1km should take 3 kwh
    max_dis = remain_kwh / 4 * 1000.0
    # max_dis = Inf
    # @show charge_time/3600.0
    # @show max_dis
    while inode <= n_node-1
        τ = τ_vec[end]
        β = β_vec[end]
        ## need to detour
        i = path[inode]
        j = path[inode+1]
        if β/B <= γ && !net.cs_flag_vec[i] && !net.cs_flag_vec[j]
            # @show inode β i j length(path_new)
            charge_time = (0.8*B-β)/0.8/1e6
            cs_idx = closest_cs(net, i, des, max_dis) 
            path_go, time_go = fastest_path(net, i, cs_idx; veh=ev, output_path=true, slow_flag=slow_flag)
            next_node_idx = closest_node(net, cs_idx, path[inode+1:end]) + inode 
            # next_node_idx = inode
            # @debug "" i cs_idx des
            next_node = path[next_node_idx]
            path_from, time_from = fastest_path(net, cs_idx, next_node ; veh=ev, output_path=true, slow_flag=slow_flag)
            path_all = vcat(path_go, path_from[2:end])
            time_all = vcat(time_go, charge_time, time_from[2:end])
            obj_tmp, β_vec_tmp, τ_vec_tmp = simulate(net, path_all, time_all, ev, β, start_time + τ2second(τ), T; check_flag=false)
            # @show path_all time_all β_vec_tmp τ_vec_tmp
            obj += obj_tmp
            t_vec_new = vcat(t_vec_new, time_all)
            path_new = vcat(path_new, path_all[1:end-1])
            β_vec = vcat(β_vec, β_vec_tmp[1:end-1])
            τ_vec = vcat(τ_vec, τ_vec_tmp[1:end-1])
            τ = τ_vec[end]
            β = β_vec[end]
            inode = next_node_idx
            i = path[inode]
            if inode == n_node
                break
            else
                j =  path[inode+1]
            end
            # continue
        end
        tij = t_vec[inode]
        road_type = get_road_type(net, i, j)
        next_b = get_next_b(net, i, j, tij, β, ev)
        next_τ = tij + τ
        push!(β_vec, next_b)
        push!(τ_vec, next_τ)
        push!(t_vec_new, tij)
        push!(path_new, i)
        inode += 1
    end
    push!(path_new, path[end])
    return path_new, t_vec_new, β_vec, τ_vec
    
end


