"""

"""

function append_graph(cs_net0::CsNet)
    n_node = length(keys(cs_net0.data))
    #cs_net = CsNet(n_node)
    #cs_net.data = deepcopy(cs_net0.data)
    cs_net = deepcopy(cs_net0)
    node_vec = Int[]
    for cs_idx in keys(cs_net.data)
        push!(node_vec, cs_idx)
        for data in cs_net.data[cs_idx]
            nei_idx = data[1]
            push!(node_vec, nei_idx)
        end
    end
    unique!(node_vec)
    for (i_id,id) in enumerate(node_vec)
        cs_net.idx2id[i_id] = id
        cs_net.id2idx[id] = i_id
        add_vertex!(cs_net.g)
    end
    # for (i_id,id) in enumerate(keys(cs_net.data))
    for id in keys(cs_net.data)
        i1 = cs_net.id2idx[id]
        #nei_ids = outneighbors(cs_net, id)
        nei_ids = [data[1] for data in cs_net.data[id]]
        for nei in nei_ids
            v = cs_net.id2idx[nei]
            Graphs.add_edge!(cs_net.g, i1, v)
        end
    end
    return cs_net
end

"""
get one edge in the cs_net. 
return: CsEdge
"""
function cs_net_get_one_edge(net::Network, src::Int, des::Int, ev::EV, dis::Real)
    # @debug "getting one edge..."
    #@timeit ep.g_to "min_cost paso" 

    # @debug "get_one_edge for src=$src, des=$des"

    if haskey(net.meta, :paso_dist_cache)
        cache = net.paso_dist_cache
    else
        cache = nothing
    end
    n_sol = 10
    step_vec = multi_paso(net, src, des, ev, ObjEnergy(), n_sol, cache)
    t_cost_vec = [(s.summary.duration, s.summary.cost) for s in step_vec]
    t_vec = [x[1] for x in t_cost_vec] 
    c_vec = [x[2] for x in t_cost_vec] 
    perm = sortperm(t_vec)
    t_vec = t_vec[perm]
    c_vec = c_vec[perm]

    min_t = minimum(t_vec)
    max_t = maximum(t_vec)
    max_cost = maximum(c_vec)
    min_cost = minimum(c_vec)
    @assert (min_t ≈ t_vec[1]) && (max_t ≈ t_vec[end]) 
    @assert (min_cost ≈ c_vec[end]) && (max_cost ≈ c_vec[1]) 
    # edge = CsEdge(min_cost, max_cost, min_t, max_t, src, des, dis, t_cost_vec)
    edge = CsEdge(min_cost, max_cost, min_t, max_t, src, des, dis, t_vec, c_vec)
    return edge
    ##################
    # step0 = paso(net, src, des, Inf; veh=ev, visitonce=false, debug_msg=false, no_charge=true)
    # min_cost = step0.summary.cost
    # max_t = step0.summary.duration
    # if min_cost > ev.cap
    #     # @debug "min_cost=$min_cost"
    #     return nothing
    # end
    # step2 = try 
    #     step2 = paso(net, src, des, 0.0; veh=ev, visitonce=false, debug_msg=false, minimize_t=true, no_charge=true, B = ev.cap)
    # catch e
    #     @error "error happens for src=$src, des=$des."
    #     rethrow(e)
    # end


    # max_cost = step2.summary.cost
    # min_t = step2.summary.duration

    # @assert (min_cost <= max_cost) 
    # @assert (min_t <= max_t) 
    # t_step = 60.0 * 30.0
    # t_vec = getrange(min_t, max_t, t_step)
    # t_cost_vec = [(min_t, max_cost)]
    # for T in t_vec[2:end-1]
    # # for cost in cost_vec
    #     # step_tmp = paso(net, src, des, 0.0; veh=ev, visitonce=true, debug_msg=false, minimize_t=true, no_charge=true, B = cost)
    #     step_tmp = paso(net, src, des, T; veh=ev, visitonce=false, debug_msg=false,  no_charge=true)
    #     tcost = (step_tmp.summary.duration, step_tmp.summary.cost)
    #     push!(t_cost_vec, tcost)
    # end
    # push!(t_cost_vec, (max_t, min_cost))
    # unique!(t_cost_vec)
    ##########
    
   
end

"""
cs_nei_flag: if true, then we need to add an edge only if there is a path not going through other charging stations.
max_neigh: Not sure when it will be useful...
"""
function compute_cs_net(net::Network, ev::EV; ow_flag = false, max_neigh::Int = typemax(Int), cs_nei_flag::Bool = false)
    cs_vec = net.cs_vec
    n = length(cs_vec)
    cs_net = load_cs_net(net.region, cs_nei_flag; tmp=true, ow_flag=ow_flag)
    
    @assert(cs_net.cs_nei_flag == cs_nei_flag)
    # min_dis = Inf
    pts = [SA[cs.lon, cs.lat] for cs in cs_vec]
    metric = Distances.Haversine()
    tree = NN.BallTree(pts, metric)
    # tree = NN.KDTree(pts)
    # node = cs_vec[1]
    cs_max_dis = j2kwh(ev.cap) / 1.5 * 1e3
    cs_min_dis = 100e3 # the minimum distance is 100 km to help reduce the graph size.
    # if net.region in ["testsmall", "test4node", "mideast"]
    # if  !(net.region in ["oldmap", "eu", ])
    if  !(net.region in ["oldmap", "eu",] || occursin("eu", net.region))
        @info "For small network, setting max/min distance to unlimited."
        cs_max_dis = Inf
        cs_min_dis = 0.0
    end
    @info "cs_max_dis=$cs_max_dis, cs_min_dis=$cs_min_dis"

    cnt = 0
    # lk = ReentrantLock()
    # max_task = Threads.nthreads()

    # if find the nearest neighbor, then augment the vector with junctions.
    if cs_nei_flag
        enum_vec = vcat(cs_vec, [ChargeStation(;idx=idx) for idx in net.junctions])
    else
        enum_vec = cs_vec
    end

    timer = time()
    t_begin = time()

    cache = ep.get_paso_dist_cache(net, ev, ObjEnergy(), 20)
    net.paso_dist_cache = cache

    # arg_vec = []
    n_computed_cs::Int = 0
    for (ics,cs) in enumerate(enum_vec)
        t_passed = time() - t_begin
        time_per_iter = t_passed / n_computed_cs
        
        @debug "ics=$ics, cs_idx=$(cs.idx)" t_passed time_per_iter n_computed_cs
        if !cs_nei_flag
            idxs, dists = NN.knn(tree, SA[cs.lon, cs.lat], round(Int,n), true)
            idx_vec_flag = ((dists .<= cs_max_dis) .&& (dists .>= cs_min_dis))
            @debug "" ics n sum(idx_vec_flag) cnt
            i_nei_vec = length(dists):-1:2
        else
            idxs = find_cs_neigh(net, cs.idx, ev)
            i_nei_vec = 1:length(idxs)
            @debug "" ics length(enum_vec) length(idxs) cnt
        end

        src = cs.idx

        arg_vec = Tuple{Int, Int, Float64}[]
        for i_nei in i_nei_vec
            idx = idxs[i_nei]
            if !cs_nei_flag
                dis = dists[i_nei]
                # if dis >= cs_max_dis || dis <= cs_min_dis || cnt >= max_neigh
                if dis >= cs_max_dis || dis <= cs_min_dis
                    continue
                end
                des = cs_vec[idx].idx
            else
                des = idx
                dis = distance3d(net, src, des)
            end
            if has_edge_net(cs_net, src, des)
                continue
            end
            # task_tmp = Threads.@spawn cs_net_get_one_edge(net, src, des, ev, dis)
            push!(arg_vec, (src, des, dis))
            # task_tmp = @task cs_net_get_one_edge(net, src, des, ev, dis)
            # push!(task_list, task_tmp)
            cnt += 1
        end #end of enumerate for the neighbors.
        if isempty(arg_vec)
            @debug "Got empty arg_vec, skip it."
            continue
        end
        
        if !haskey(cs_net.data, src)
            n_computed_cs += 1
        else
            t_passed = 0.0
        end
        cs_net = _compute_one_cs!(net, ev, cs_net, arg_vec)

        time_limit = 300 # 5 min to save the results. To reduce saving overhead.
        if (time() - timer) > time_limit
            timer = time()
            save_cs_net(net.region, cs_net, cs_nei_flag; tmp=true)
            GC.gc(false)
        end
    end
    save_cs_net(net.region, cs_net, cs_nei_flag; tmp=true)
    cs_net1 = append_graph(cs_net)
    return cs_net1
end

"""
compute the edges for one charging station
arg_vec: a vector of (src, des, dis) tuple
"""
function _compute_one_cs!(net::Network, ev::EV, cs_net::CsNet, arg_vec::Vector{T}) where T <: Any
    # gradually schedule and yield the task 
    ntask = length(arg_vec)
    GC.gc(false)

    # nthread = Threads.nthreads()
    # edge_list = Vector(undef, ntask)
    # edge_list = []
    # pbar = ProgressMeter.Progress(ntask; dt=0.01, showspeed=true)
    @debug "computing for ntask=$ntask"
    function func(arg)
        (src, des, dis) = arg
        return cs_net_get_one_edge(net, src, des, ev, dis) 
    end
    # thunk_size = 
    edge_list = ep.e2map(func, arg_vec, true, false; )
    # task_list = [@task cs_net_get_one_edge(net, src, des, ev, dis)  for (src, des, dis) in arg_vec]
    #Threads.@threads for arg in arg_vec
    #    (src, des, dis) = arg
    #    edge = cs_net_get_one_edge(net, src, des, ev, dis) 
    #    lock(lk) do 
    #        if !isnothing(edge)
    #            push!(edge_list, edge)
    #        end
    #        ProgressMeter.next!(pbar)
    #        # flush(pbar.output)
    #    end
    #end
    #ProgressMeter.finish!(pbar)
    # ntask = length(task_list)
    # ntask_each_thread = zeros(nthread)
    # task_thread_id = zeros(Int, ntask)
    # ntask_each_thread[1] = Inf
    # for (itask,t) in enumerate(task_list)
    #     # t may now get run on any thread
    #     # task.stikcy means we statically set its 
    #     t.sticky = true
    #     # ntask_done = sum(istaskdone.(task_list))
    #     # ntask_run = itask - ntask_done
    #     # update the number of task running for each thread
    #     for i in 2:nthread
    #         ntask_each_thread[i] = 0
    #         for id in task_thread_id
    #             if id != 0
    #                 ntask_each_thread[id] += 1
    #             end
    #         end
    #         # The follow is not true, an scheduled task can wait in the queue and not get started.
    #         # occur_vec= findall(x->(Threads.threadid(x)==i && !istaskdone(x) && istaskstarted(x)), task_list)
    #     end
    #     # @show ntask_each_thread

    #     # if the number of running task is larger than cpu cores, yield the current core and wait for available cores.
    #     if all(ntask_each_thread .>= 4) || nthread == 1
    #         cur_id = Threads.threadid()
    #         set_task_tid(t, cur_id)
    #         schedule(t)
    #         # @debug "yielding with" 
    #         # @show ntask_each_thread
    #         task_thread_id[itask] = cur_id
    #         yield()
    #         # yieldto(t)
    #     else
    #         # find the thread with minimum number of running task.
    #         # @debug "scheduling with " 
    #         # @show ntask_each_thread
    #         idx = argmin(ntask_each_thread)
    #         # @show itask idx
    #         task_thread_id[itask] = idx
    #         set_task_tid(t, idx)
    #         schedule(t)
    #     end
    #     GC.safepoint()

    #     # if task is done, we set its thread id to zero.
    #     task_thread_id[istaskdone.(task_list)] .= 0
    #     ntask_done = sum(istaskdone.(task_list))
    #     # ntask_fail = sum(istaskfailed.(task_list))
    #     # ntask_start = sum(istaskstarted.(task_list))
    #     # @show ntask_done,ntask_fail,ntask_start
    #     ProgressMeter.update!(pbar, ntask_done)
    #     # force to print the current progress
    #     flush(pbar.output)
    # end
    # edge_list = fetch.(task_list)
    # ProgressMeter.finish!(pbar)

    @debug "Got $(length(edge_list)) edges."
    if !isempty(edge_list)
        for edge in edge_list
            if !isnothing(edge)
                add_edge_net!(cs_net, edge.src, edge.des, edge)
            end
        end
    else
        @warn "get empty edge_list" 
    end
    # @debug "Sleep 1 seconds for interrupt..."
    # sleep(1)
    return cs_net
end

"""
src: the index of a cs
"""
function find_cs_neigh(net::Network, src::Int, ev::EV)
    n_node = nv(net)
    q = PriorityQueue{Int, Float64}()
    (;cs_vec) = net
    cs_idx_vec = [cs.idx for cs in cs_vec]
    cs_idx_r_vec = [cs.idx_r for cs in cs_vec]
    # cs_idx_vec = vcat(cs_idx_vec, net.junctions, cs_idx_r_vec)
    # bfs to find the neighbor of a charging station.
    q[src] = 0.0
    
    visited = falses(n_node)
    visited[src] = true

    cs_flag = falses(n_node)
    cs_flag[cs_idx_vec] .= true
    cs_idx_r_flag = falses(n_node)
    cs_idx_r_flag[cs_idx_r_vec] .= true
    junc_flag = falses(n_node)
    junc_flag[net.junctions] .= true

    g_score = fill(Inf, n_node)
    g_score[src] = 0

    nei_vec = Int[]
    while !(isempty(q)) 
        u = dequeue!(q)
        visited[u] = true
        @inbounds for v in outneighbors(net, u)
            if !visited[v]
                # @debug "visiting $u -> $v"
                #min_t, max_t = minmax_t(net, u, v, ev; nowait=true)
                min_spd, max_spd = get_minmax_speed(net, u, v, ev)
                min_fuelcost = cost(ObjEnergy(), ev, net, u, v, min_spd)
                g_score_tmp = g_score[u] + min_fuelcost
                if g_score_tmp > ev.cap
                    continue
                end
                if (cs_flag[v] || junc_flag[v])
                    push!(nei_vec, v)
                    continue
                end
                # if v is the real world node
                if cs_idx_r_flag[u] && (v < net.cs_start_idx)
                    continue 
                end
                if g_score_tmp < g_score[v]
                    g_score[v] = g_score_tmp
                    q[v] = g_score_tmp
                    # if cs_flag[v]  
                    #     ## NOTE: this will reudce the edge for the case where two charging stations are neighbors in the original graph.
                    #     if !cs_idx_r_flag[v]
                    #         push!(nei_vec, v)
                    #     end
                    # else
                    #     q[v] = g_score_tmp
                    # end
                end
            end
        end
    end
    unique!(nei_vec)
    return nei_vec
end

"""
modify the max and min cost to satisfy the capacity limit
"""
function cs_net_restrict_cap!(cs_net::CsNet, B::Float64)
    # for  
    for edge_p_vec in values(cs_net.data)
        # edge_p = (des, edge)
        for (i_edge_p, edge_p) in enumerate(edge_p_vec)
            edge = edge_p[2]
            # (;t_cost_vec) = edge
            (;t_vec, c_vec) = edge
            # cost_vec = [x[2] for x in t_cost_vec]
            # @show edge.src edge.des cost_vec

            if isempty(c_vec)
                continue
            end

            if c_vec[end] > B
                edge_p_vec[i_edge_p] = (edge_p[1], CsEdge(;src=edge.src,des=edge.des))
                continue
            end


            idx::Int = length(t_vec)
            for idx0 in 1:length(t_vec)
                if c_vec[idx0] <= B
                    idx = idx0
                    break
                end
            end

            # t_cost_vec1 = t_cost_vec[1:idx]
            t_vec1::Vector{Float64} = t_vec[idx:end]
            c_vec1::Vector{Float64} = c_vec[idx:end]
            @assert(maximum(c_vec1) <= B)

            # t_vec1 = [x[1] for x in t_cost_vec1]
            # c_vec1 = [x[2] for x in t_cost_vec1]
            min_t = minimum(t_vec1)
            max_t = maximum(t_vec1)
            min_cost = minimum(c_vec1)
            max_cost = maximum(c_vec1)
            # searchsortedafter
            new_edge = CsEdge(min_cost, max_cost, min_t, max_t, edge.src, edge.des, edge.dis, t_vec1, c_vec1)
            # edge_p[2] = new_edge
            edge_p_vec[i_edge_p] = (edge_p[1], new_edge)
        end
    end
    return cs_net
end