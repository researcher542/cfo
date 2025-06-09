

"""
n1, n2, n_mid are nodes in the network.

The angle for <n1, n_mid, n2> ∈ [0, π].
"""
function AngleBetweenVectors.angle(n1, n2, n_mid)
    v1 = (n1.lat - n_mid.lat, n1.lon - n_mid.lon)
    v2 = (n2.lat - n_mid.lat, n2.lon - n_mid.lon)
    return angle(v1, v2)
end

"""
A heuristic to find the "best" charging station.
"""
function find_charge_station(net::Network, src::Int, des::Int, β0::Real, kwh_per_km)
    src_node = getnode(net, src)
    des_node = getnode(net, des)
    # Based on our poster, the e-truck uses 2 kwh/km for t-ratio=1.05
    # 
    max_dis = j2kwh(β0) / kwh_per_km * 1e3
    cs_can_vec = []
    for cs in net.cs_vec
        node = getnode(net, cs.idx)
        dis = distance3d(node, src_node)
        dis2des = distance3d(node, des_node)
        # The angle, close to zero is better.
        # ang = angle(src_node, des_node, node) / pi * 180
        # degree of 30
        # ang_flag = ang <= 30
        dis_flag = (dis <= max_dis)  && (dis >= 0.6 * max_dis)
        if dis_flag
            # @debug "" ang_flag dis_flag ang dis
        end
        score = -dis + dis2des
        price = g_mean_price_dict[cs.region]
        if dis_flag
            push!(cs_can_vec, (cs, dis, price, score))
        end
    end
    if isempty(cs_can_vec)
        @exfiltrate
        @assert false
    end
    # price_v = (data[3] for data in cs_can_vec)
    # min_price = minimum(price_v)
    # idx = findall(x->x[3]==min_price, cs_can_vec)
    # cs_can_vec = cs_can_vec[idx]
    val, idx = findmin([data[4] for data in cs_can_vec])
    return cs_can_vec[idx][1]
end

"""
A greedy heuristic to get a feasible solution.
"""
function cfo_heuristic(net::Network, src::Int, des::Int, ev::EV, start_time::DateTime, T::Real, β0::Real)
    kwh_per_km = 4
    max_dis = j2kwh(β0) / kwh_per_km * 1e3
    cur_node_idx = src
    dis = distance3d(net, cur_node_idx, des)
    cur_β = β0
    path = Vector{Int}([src])
    t_vec = Vector{Float64}()
    wait_idx_vec = Vector{Int}()
    mid_cs_vec = Vector{ChargeStation}()
    B = ev.cap
    cnt = 0
    obj = 0.0
    while dis > max_dis
        @debug "" cnt dis/1e3 max_dis/1e3 cur_node_idx 
        cs = find_charge_station(net, cur_node_idx, des, cur_β, kwh_per_km)
        # path_go, time_go = fastest_path(net, cur_node_idx, cs.idx; veh=ev, output_path=true, slow_flag=)
        path_go, time_go = paso(net, cur_node_idx, cs.idx, 4*3600.0; veh=ev, output_step=false, visitonce=true)
        next_node_idx = outneighbors(net, path_go[end])[1]

        τ = sum(t_vec)
        obj_tmp, β_vec_tmp, τ_vec_tmp = simulate(net, path_go, time_go, ev, cur_β, start_time + τ2second(τ), T; check_flag=false)
        obj += obj_tmp
        # only charge to 80 %
        cur_β = β_vec_tmp[end]
        charge_time = (0.8*B-cur_β)/1e6
        path = vcat(path, path_go[2:end], next_node_idx)
        t_vec = vcat(t_vec, time_go, charge_time)

        push!(mid_cs_vec, cs)
        push!(wait_idx_vec, length(path)-2)
        cur_β = get_next_b(net, path_go[end], next_node_idx, charge_time, cur_β, ev)
        max_dis = j2kwh(cur_β) / kwh_per_km * 1e3
        cur_node_idx = next_node_idx
        dis = distance3d(net, cur_node_idx, des)
        cnt += 1
    end
    remain_time = T - sum(t_vec)


    path_go, time_go = fastest_path(net, cur_node_idx, des; veh=ev, output_path=true, slow_flag=false)
    # next_node_idx = outneighbors(net, path_go[end])[1]
    path = vcat(path, path_go[2:end])
    t_vec = vcat(t_vec, time_go)

    return path, t_vec
end
