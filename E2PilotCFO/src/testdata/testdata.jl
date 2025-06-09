


"""
A high-level function to store some predefined test problems.
"""
function get_test_prob(region::String; kwargs...)
    if region == "ATNA"
        prob = get_test_prob_us_reduced(region; box_tol=0.0, merge_dis=5_000.0, kwargs...)
    elseif region == "LASD"
        prob = get_test_prob_us_reduced(region; box_tol = 0.0, merge_dis=5_000.0, B = 200.0 * 3.6e6, kwargs...)
    elseif region == "TXOH"
        prob = get_test_prob_us_reduced(region; box_tol=1.0, kwargs...)
    elseif region == "oldmap"
        prob = get_test_prob_large(;kwargs...)
    elseif region == "eu"
        prob = get_test_prob_large("eu.merged.v2"; 
            carbon_dataset = ElecMapDataset(), 
            sparse_cs_dis = [5_000.0, 20_000.0][2], 
            kwargs...)
    else
        @error "Unknown region $region"
        return missing
    end
    return prob
end

function get_test_src_des_latlon(region, verbose=false)
    od_str_dict = OrderedDict(
        "TXOH" =>  ( "484-Dallas-Fort Worth TX-OK (TX Part)", "391-Cincinnati OH-KY-IN (OH Part)"),
        "ATNA" => ("131-Atlanta GA", "472-Nashville TN"),
        "LASD" => ("061-Los Angeles CA", "063-San Diego CA"),
        "test" => ("061-Los Angeles CA", "063-San Diego CA"),

        ### The following are selected values
        "LACL" => ("061-Los Angeles CA", ""),
    )
    if haskey(od_str_dict, region)
        (src_str, des_str) = od_str_dict[region]
    else
        ## region should have the format of "srcidx-desidx"
        parts = split(region, "-")
        src_idx = parse(Int, parts[1])
        des_idx = parse(Int, parts[2])
        src_str = ep.idx2fafregion(src_idx) 
        des_str = ep.idx2fafregion(des_idx)
        if src_str == "" || des_str == ""
            @error "Unknown region $region"
            return missing, missing
        end
        if verbose 
            @show src_idx, des_idx
        end
    end
    if verbose
        @show src_str, des_str
    end

    faf_reg_dict = ep.read_faf_region()
    src_latlon = faf_reg_dict[src_str]
    des_latlon = faf_reg_dict[des_str]
    if verbose 
        @show ep.distance2d(src_latlon, des_latlon) / ep.g_meter_per_mile
    end
    return src_latlon, des_latlon
end


"""
Get the synthetic test problem with 4 nodes.
"""
function get_test_prob_small(
    region::String;
    st = g_default_st,
    carbon_dataset = CambiumDataset(),
    T = 2*3600.0,
    scenario = g_default_scenario,
    sep_cs_node_flag = false,
    cs_nei_flag = false,
    predict_mode = PredictPerfect(),
    B = 3.6e6 * 100.0,
    )

    net = ep.readmapdata(region; use_saved=false)
    src,des = 1,4
    net.sep_cs_node_flag = sep_cs_node_flag
    cs_vec = ChargeStation[]
    for idx in 1:nv(net)
        node = ep.getnode(net, idx)
        cs = ChargeStation(lat=node.lat, lon=node.lon, idx=idx, region="CA")
        push!(cs_vec, cs)
    end
    ev = get_ev()


    add_charge_station!(net, cs_vec)

    @timeit ep.g_to "load_cs_net" cs_net = load_cs_net(net.region, cs_nei_flag)
    ev.cap = B

    cs_net = cs_net_restrict_cap!(cs_net, ev.cap)
    net.cs_net = cs_net
    carbon_dict = get_carbon_dict(carbon_dataset)
    carbon_predict_dict = get_carbon_dict_predict(carbon_dataset)
    
    fix_data = CfoProbData(;
        net=net, ev=ev, 
        cs_dict=net.cs_dict, 
        carbon_dict=carbon_dict, 
        carbon_predict_dict=carbon_predict_dict)
    src,des = 1, 4

    prob = CfoProb(;
        fix_data=Ref(fix_data),
        src=src, des=des, 
        β0=ev.cap, T=T, 
        start_time=st, 
        scenario=scenario,
        predict_mode=predict_mode,
        odset = FAF()
    )
    return prob
    
end

function get_latlon_box(net, nd_vec)
    min_lat = minimum( ep.getnode(net, i).lat for i in nd_vec )
    max_lat = maximum( ep.getnode(net, i).lat for i in nd_vec )
    min_lon = minimum( ep.getnode(net, i).lon for i in nd_vec )
    max_lon = maximum( ep.getnode(net, i).lon for i in nd_vec )
    return min_lat, max_lat, min_lon, max_lon
end

function get_reduced_box(net0::AbsNet, src0, des0, box_tol, veh)
    ## Make share the shortest from src to des is within the reduced network.
    # path0 = ep.shortest_path(ep.Astar(), net0, src0, des0; throw_err_flag=false)

    # distmx = zeros(nv(net0), nv(net0))
    # I = Int[]
    # J = Int[]
    # V = Float64[]
    # for (k, edge) in net0.waydata
    #     # dis = edge.distance
    #     # distmx[edge.src, edge.des] = dis
    #     push!(I, edge.src)
    #     push!(J, edge.des)
    #     push!(V, edge.distance)
    # end
    # @time distmx = sparse(I, J, V)

    path0 = ep.shortest_path(ep.Astar(), net0, src0, des0; throw_err_flag=false)
    path1 = ep.shortest_path(ep.Astar(), net0, src0, des0; getdist = (net, u, v, lambda) -> minmax_t(net, u, v, veh)[1] ,)
    path2 = ep.shortest_path(ep.Astar(), net0, src0, des0; getdist = (net, u, v, lambda) -> minmax_t(net, u, v, veh)[2] ,)

    all_nd_vec = vcat(path0, path1, path2)
    # all_nd_vec = vcat(path0, )

    # K = 2
    # @time yen_state = Graphs.yen_k_shortest_paths(net0.g, src0, des0, distmx, K) 
    # path_vec = yen_state.paths
    # all_nd_vec = vcat(path_vec...)
    (min_lat, max_lat, min_lon, max_lon) = get_latlon_box(net0, all_nd_vec)

    red_net_kwargs = (
            min_lat= min_lat - box_tol, 
            max_lat= max_lat + box_tol, 
            min_lon= min_lon - box_tol, 
            max_lon= max_lon + box_tol, 
            )
    
    return red_net_kwargs
end


function reduced_idxmap(net0::AbsNet, net1::AbsNet, idx_vec::Vector{Int})
    return [reduced_idxmap(net0, net1, idx) for idx in idx_vec]
end

"""
Convert the node from net0 to net1
"""
function reduced_idxmap(net0::AbsNet, net1::AbsNet, idx::Int)
    nd0 = ep.getnode(net0, idx)
    nd1_idx = ep.closest_node_in_max_com(net1, LatLon(nd0.lat, nd0.lon))

    nd1 = ep.getnode(net1, nd1_idx)
    @assert nd0.lat == nd1.lat && nd0.lon == nd1.lon (@exfiltrate; "The node $(nd0) in the original network is not the same as the node $(nd1) in the reduced network $(net1.region).")

    return nd1_idx 
end

"""
Reduce the large cs_net to the reduced network.

cs_vec: the charge stations in the reduced network.
"""
function reduce_cs_net(net0::AbsNet, net1::AbsNet, cs_vec, cs_net0)
    cs_net1 = CsNet(-1; cs_net0.cs_nei_flag)

    # in_reduced_net = falses(nv(net0))
    # for nd in net1.nodesdata
    #     idx = ep.closest_node(net0, LatLon(nd.lat, nd.lon))
    #     in_reduced_net[idx] = true
    # end
    # in_cs_net = falses(nv(net0))
    idx_dict = Dict{Int, Int}()
    for cs in cs_vec
        idx = ep.closest_node_in_max_com(net0, LatLon(cs.lat, cs.lon))
        nd = ep.getnode(net0, idx)
        @assert nd.lat == cs.lat && nd.lon == cs.lon "The node $(nd) in the original network is not the same as the charge station $(cs) in the reduced network $(net1.region)."
        idx_dict[idx] = cs.idx
    end

    for (o, edge_vec) in cs_net0.data
        if !haskey(idx_dict, o)
            continue
        end
        for (d, edge) in edge_vec
            if !haskey(idx_dict, d)
                continue
            end
            # @show edge.src edge.des
            # @show idx_dict[edge.src] idx_dict[edge.des]
            new_src = reduced_idxmap(net0, net1, edge.src)
            new_des = reduced_idxmap(net0, net1, edge.des)
            new_edge = CsEdge(
                src = new_src,
                des = new_des,
                min_cost = edge.min_cost,
                max_cost = edge.max_cost,
                min_t = edge.min_t,
                max_t = edge.max_t,
                dis = edge.dis,
                t_vec = edge.t_vec,
                c_vec = edge.c_vec,
            )
            add_edge_net!(cs_net1, new_edge.src, new_edge.des, new_edge)
        end
    end

    # for (ics, cs) in enumerate(cs_vec)
    #     add_edge_net!(cs_net1, edge.src, edge.des, edge)
    # end
    cs_net2 = append_graph(cs_net1)
    return cs_net2
end

"""
Get a test problem from the U.S. dataset with reduced nodes and edges.
"""
function get_test_prob_us_reduced(
    region::String,
    ;
    box_tol::Float64 = 1e-1, ## 1 degree is about 100 km
    st = g_default_st,
    scenario = g_default_scenario,
    predict_mode = PredictPerfect(),
    carbon_dataset = CambiumDataset(),
    sep_cs_node_flag = false,
    cs_nei_flag = false,
    T = 0.0,
    B = 3.6e6 * 1000.0,
    merge_dis = 0.0,
    test_od_flag::Bool = true, ## Test if OD is connected in the network
    sparse_cs_dis::Real = 5_000.0,
    )
    net0 = ep.readmapdata("oldmap"; use_saved=false)
    src_latlon, des_latlon = get_test_src_des_latlon(region)

    src0 = ep.closest_node_in_max_com(net0, src_latlon)
    des0 = ep.closest_node_in_max_com(net0, des_latlon)


    path0 = ep.shortest_path(ep.Astar(), net0, src0, des0; throw_err_flag=false)

    ev = get_ev()

    red_net_kwargs = get_reduced_box(net0, src0, des0, box_tol, ev)
    
    @show src_latlon, des_latlon red_net_kwargs
    @info "reducing to $region"
    net_red1 = ep.reduce_net(net0, region; red_net_kwargs...)

    if merge_dis > 0.0
        @info "merging nodes with distance $merge_dis"
        net_red1 = ep.net_merge_nodes(net_red1, region, merge_dis) 
    end 
    @info "For reduced net: nv=$(nv(net_red1)) ne=$(ne(net_red1))"

    net = net_red1


    src = ep.closest_node_in_max_com(net, src_latlon)
    des = ep.closest_node_in_max_com(net, des_latlon)

    cs_vec0 = get_charge_station(net0, scenario, fuel=false)

    ### Note that previously sparse_cs_dis is 20_000.0. We make it to 5000.0 to make it consistent with the oldmap dataset.
    cs_vec0 = sparsify_cs_vec(net0, cs_vec0, sparse_cs_dis)
    cs_vec = Vector{ChargeStation}()
    for cs in cs_vec0
        node_idx = ep.closest_node_in_max_com(net, LatLon(cs.lat, cs.lon))
        cs_new = deepcopy(cs)
        nd_red = getnode(net, node_idx)
        if nd_red.lat == cs.lat && nd_red.lon == cs.lon
            cs_new.idx = node_idx
            push!(cs_vec, cs_new)
        end
    end

    net.sep_cs_node_flag = sep_cs_node_flag
    add_charge_station!(net, cs_vec)
    cs_net = load_cs_net(net0.region, cs_nei_flag)
    ev.cap = B
    cs_net = cs_net_restrict_cap!(cs_net, ev.cap)

    cs_net1 = reduce_cs_net(net0, net, cs_vec, cs_net)

    net.continent = :us
    net.cs_net = cs_net1
    beta0 = ev.cap
    @show length(cs_vec)
    @show nv(net)
    # energy cost [1813kWh, 2.088kWh], distance ~1179km time [13.9h, 24.5h]

    src = ep.closest_node_in_max_com(net, LatLon(src_latlon.lat + 0.0, src_latlon.lon + 0.0) )
    des = ep.closest_node_in_max_com(net, LatLon(des_latlon.lat + 0.0, des_latlon.lon - 0.0) )

    ## check if the reduced network remain the same shortest path as the original network.
    path1 = ep.shortest_path(ep.Astar(), net, src, des)
    for i in 1:length(path1)-1
        nd0 = ep.getnode(net0, path0[i])
        nd1 = ep.getnode(net, path1[i])
        # @show nd0 nd1
        @assert nd0.lat == nd1.lat && nd0.lon == nd1.lon "The node $(i) in the path is not the same: $(nd0) vs $(nd1)"
    end

    ## another set of check
    res_file_path = joinpath(g_cfo_result_dir, "idx0.fast-ds.$(region).jld2")
    if isfile(res_file_path) && false
        res = load_object(res_file_path)
        for sol in res.primal.sub_sol_vec
            for ind in 1:length(sol.path)-1
                nd01 = ep.getnode(net0, sol.path[ind])
                nd02 = ep.getnode(net0, sol.path[ind+1])

                edge0 = ep.getway(net0, nd01.idx, nd02.idx)

                nd1_idx = ep.closest_node_in_max_com(net, LatLon(nd01.lat, nd01.lon))
                nd2_idx = ep.closest_node_in_max_com(net, LatLon(nd02.lat, nd02.lon))

                nd11 = ep.getnode(net, nd1_idx)
                nd12 = ep.getnode(net, nd2_idx)
                
                if nd11.idx == nd12.idx
                    continue
                end
                edge1 = ep.getway(net, nd11.idx, nd12.idx)

                @assert nd01.lat == nd11.lat && nd01.lon == nd11.lon "The node $(nd01) in the path is not the same: $(nd01) vs $(nd11)"
                @assert edge0.distance == edge1.distance
                # @assert nd0.lat == nd1.lat && nd0.lon == nd1.lon "The node $(nd) in the path is not the same: $(nd0) vs $(nd1)"
            end
        end
    end
    
    ## shift the source node to north to let the default route has high carbon intensity.
    # des = ep.closest_node_in_max_com(net, LatLon(des_latlon.lat + 0.0, des_latlon.lon - 1.0) )

    src = closest_cs_node(net, src)
    des = closest_cs_node(net, des)

    carbon_dict = get_carbon_dict(carbon_dataset)
    carbon_dict_predict = get_carbon_dict_predict(carbon_dataset)
    fix_data = CfoProbData(;
        net=net, ev=ev, 
        cs_dict=net.cs_dict, 
        carbon_dict=carbon_dict, 
        carbon_predict_dict=carbon_dict_predict)

    

    path = ep.shortest_path(ep.Astar(), net, src, des; throw_err_flag=false)

    if isempty(path)
        @warn "The src and des are not connected in the network. We might need to increase the box_tol."
    else
        if test_od_flag
            @info "The src and des are connected in the network." length(path)
        end
    end

    if T == 0.0
        total_dis = sum(distance3d(net, path[i], path[i+1]) for i in 1:length(path)-1)
        speed = 10.0
        T = total_dis / speed
        ## Add the consideration of charing time
        T += total_dis / (300 * 1e3) * 3600.0
        @debug "Setting the deadline to $(T/3600.0) hours."
    end
    
    prob = CfoProb(;
        fix_data=Ref(fix_data),
        src=src, des=des, β0=beta0, T=T, 
        start_time=st, 
        scenario=scenario,
        predict_mode=predict_mode,
        odset = FAF(),
    )
    
    return prob
end

"""
Check if a problem dataset is valid.
Check if the dataset is consistent for regions, od pairs, etc.
"""
function check_prob(prob)
    continent = prob.net.continent
    region = prob.net.region
    (; predict_mode, odset, carbon_dataset) = prob

    @assert ep.region2continent(region) == continent

    if continent == :us
        @assert odset == FAF()
    elseif continent == :eu
        @assert predict_mode == PredictPerfect()
        @assert odset == ETIS()
        @assert carbon_dataset == ElecMapDataset() 
    end
end

"""
The major U.S. data we used in the paper. Contains about 80,000+ nodes.
"""
function get_test_prob_large(
    region::String = "oldmap",
    ;
    carbon_dataset = CambiumDataset(),
    scenario = g_default_scenario,
    sep_cs_node_flag = false,
    cs_nei_flag = false,
    predict_mode = PredictPerfect(),
    T::Real = 60*3600.0,
    src::Int = -1,
    des::Int = -1,
    sparse_cs_dis::Real = 5_000.0,
    )

    continent = ep.region2continent(region)
    @show carbon_dataset
    if continent == :eu
        odset = ETIS()
        if src == -1 || des == -1
            src, des = (121727, 79859) 
        end
    elseif continent == :us
        odset = FAF()
        if src == -1 || des == -1
            src, des = (17745, 81148) 
        end
    else
        throw(ArgumentError("Unknown region $region"))
    end
    net = ep.readmapdata(region; use_saved=false)
    @timeit g_to "get_cs_vec" cs_vec = get_charge_station(net, scenario; fuel=false, carbon_dataset=carbon_dataset, continent = continent)

    cs_vec = sparsify_cs_vec(net, cs_vec, sparse_cs_dis)

    net.sep_cs_node_flag = sep_cs_node_flag
    add_charge_station!(net, cs_vec)
    # @show length(cs_vec)

    cs_net = load_cs_net(net.region, cs_nei_flag)
    ev = get_ev()
    cs_net = cs_net_restrict_cap!(cs_net, ev.cap)
    net.cs_net = cs_net
    net.cs_vec = cs_vec

    # cs_vec_idx = [cs.idx for cs in cs_vec]
    carbon_dict = get_carbon_dict(carbon_dataset)
    carbon_dict_predict = get_carbon_dict_predict(carbon_dataset)
    start_time = g_default_st

    fix_data = CfoProbData(;
        net=net, ev=ev, 
        cs_dict=net.cs_dict, 
        carbon_dict=carbon_dict, 
        carbon_predict_dict=carbon_dict_predict,
        carbon_dataset = carbon_dataset,
    )
    
    src = closest_cs_node(net, src)
    des = closest_cs_node(net, des)
    beta0 = ev.cap

    prob = CfoProb(;
        fix_data=Ref(fix_data),
        src=src, des=des, β0=beta0, T=T, 
        start_time=start_time, 
        scenario=scenario,
        predict_mode=predict_mode,
        odset = odset,
    )

    check_prob(prob)
    return prob
end


"""
A function to generate test data for the Gurobi solver.
"""
function generate_grb_test_data(region::String = "grb.test")
    src_str, des_str = [
    ("132-Savannah GA", "242-Washington DC-VA-MD-WV (MD Part)"),
    ("131-Atlanta GA", "472-Nashville TN")
    ][2]
    faf_reg_dict = ep.read_faf_region()
    src_latlon = faf_reg_dict[src_str]
    des_latlon = faf_reg_dict[des_str]
    net0 = ep.readmapdata("oldmap"; use_saved=false)
    tol = 0.0
    kwargs = (
              min_lat=min(src_latlon.lat, des_latlon.lat) - tol, 
              max_lat=max(src_latlon.lat, des_latlon.lat) + tol, 
              min_lon=min(src_latlon.lon, des_latlon.lon) - tol, 
              max_lon=max(src_latlon.lon, des_latlon.lon) + tol, 
            )
    @info "reducing to mideast"
    net_red1 = ep.reduce_net(net0, region; kwargs...)
    @info "merging nodes"
    net_red2 = ep.net_merge_nodes(net_red1, region, 1_000.0)
    net = net_red2
    src = ep.closest_node_in_max_com(net, src_latlon)
    des = ep.closest_node_in_max_com(net, des_latlon)

    @timeit ep.g_to "get_cs_vec" cs_vec = get_charge_station(net, fuel=false)
    cs_vec = sparsify_cs_vec(net, cs_vec, 5_000.0, [src, des])
    net.sep_cs_node_flag = false
    add_charge_station!(net, cs_vec)
    @timeit ep.g_to "load_cs_net" cs_net = load_cs_net(net.region, false)
    @timeit ep.g_to "get_ev" ev = get_ev()
    ev.cap = 3.6e6 * 350
    cs_net = cs_net_restrict_cap!(cs_net, ev.cap)
    net.cs_net = cs_net

    β0 = ev.cap

    src0, des0 = src, des
    src = closest_cs_node(net, src)
    des = closest_cs_node(net, des)
    start_time = DateTime(2021, 4, 4, 3, 0, 0)
    T = 3600.0 * 24.0 * 2.0


    prob = CfoProb(;net=net, ev=ev, src=src, des=des, β0=β0, T=T, start_time=start_time, cs_dict=net.cs_dict)

    return prob
end