function _remove_invalid_candidates!(net::AbsNet, tree, node_idx::Int, candidate_set::Set{Int}, candidate_vec::Vector{Int}, min_dist::Float64, npts::Int)
    k = 5
    cur_nd = getnode(net, node_idx)
    cur_pt = SA[cur_nd.lat, cur_nd.lon]
    dists = zeros(k)
    idxs = zeros(Int, k)
        
    while dists[end] < min_dist && k < npts
        k = min(k * 2, npts)
        idxs, dists = NN.knn(tree, cur_pt, k, )
        sorted_idx = sortperm(dists)
        idxs = idxs[sorted_idx]
        dists = dists[sorted_idx]
    end
        
    for (idx, dist) in zip(idxs, dists)
        if dist < min_dist
            nd = candidate_vec[idx]
            delete!(candidate_set, nd)
        else
            break
        end
    end
end

"""
Merge two nodes if there are too close
"""
function net_merge_nodes(net::Network{NT, WT}, region::String, min_dist::Float64) where {NT <: AbstractNode, WT <: AbstractWay}

    function node_old2new(net, new_net, old_idx::Int)
        node = getnode(net, old_idx)
        id = node.id
        new_idx = get(new_net.id2idx, id, -1)
        return new_idx
    end
    # net = deepcopy(net0)
    add_junction!(net; min_nei=3)

    candidate_vec = copy(net.max_com)
    candidate_set = Set(net.max_com)
    new_node_vec = intersect(net.junctions, candidate_set)
    # new_node_vec = Int[]
    @debug "Searching for new nodes with $(length(new_node_vec)) junctions, and max_com = $(length(net.max_com))"

    # for (i, node_idx) in enumerate(net.max_com)
    # for (i, node_idx) in enumerate(candidate_vec)

    pts = [SA[getnode(net, nd).lat, getnode(net, nd).lon] for nd in candidate_vec]
    metric = Distances.Haversine()
    tree = NN.BallTree(pts, metric)
    i = 0
    npts = length(candidate_vec)
    for node_idx in new_node_vec
        _remove_invalid_candidates!(net, tree, node_idx, candidate_set, candidate_vec, min_dist, npts)
    end

    while !isempty(candidate_set)
        node_idx = pop!(candidate_set)
        i += 1
        if is_endpoint(net, node_idx) # || !(node_idx in net.max_com)
            continue
        end
        if i % 1000 == 0
            @debug "processing $(i)-th n_node=$(length(new_node_vec)) total $(length(candidate_set))"
        end
        # add_flag = true
        ### NOTE: we do not need to check the distance between the new node and the existing nodes as we have removed them in the candidate_vec
        
        _remove_invalid_candidates!(net, tree, node_idx, candidate_set, candidate_vec, min_dist, npts)
        push!(new_node_vec, node_idx)
        

        # dis_vec = zeros(length(candidate_vec))
        # n = length(candidate_vec)
        # candidate_vec1 = collect(candidate_vec)
        # # Threads.@threads 
        # @batch for i in 1:n
        #     nd = candidate_vec1[i]
        #     dis = distance3d(net, node_idx, nd)
        #     dis_vec[i] = dis
        # end
        # to_remove_candidate_vec = Int[]
        # for i in 1:n
        #     dis = dis_vec[i]
        #     if dis < min_dist
        #         push!(to_remove_candidate_vec, candidate_vec1[i])
        #     end
        # end
        # for nd in to_remove_candidate_vec
        #     delete!(candidate_vec, nd)
        # end
    end

    @debug "new node vec $(length(new_node_vec))"

    # create the new network
    new_net = Network(region)
    node_cnt = 1
    # @debug "adding new node"
    for node_idx in new_node_vec
        old_node = getnode(net, node_idx)
        (;lat,lon,ele,id) = old_node
        new_node = NT(;lat=lat, lon=lon, ele=ele, idx=node_cnt, id=id)
        addnode!(new_net, new_node)
        node_cnt += 1
    end

    way_cnt = 1

    # @debug "adding new ways"
    new_node_set = Set(new_node_vec)
    for u in new_node_vec
        # @debug "adding ways for $u"
        u_new_idx = node_old2new(net, new_net, u)
        neighs = find_neigh_bfs(net, u, new_node_set)
       

        # u_node = getnode(net, u)
        for v in neighs
            v_new_idx = node_old2new(net, new_net, v)
            if has_way(new_net, u_new_idx, v_new_idx)
                continue
            end
            # v_node = getnode(net, v)
            # theta = grade(u_node, v_node)
            # dis = distance3d(u_node, v_node)
            
            # new_way = WT(; src=u_new_idx, des=v_new_idx, id=way_cnt, nodes=[u_new_idx, v_new_idx], grades=[theta], distances=[dis])
            new_way = WT(net, u_new_idx, v_new_idx, way_cnt, false, [u_new_idx, v_new_idx], )
            addway!(new_net, new_way)
            uu = inneighbors(net, v)[1]
            spd = net.speeddata[(uu, v)]
            new_net.speeddata[(new_way.src, new_way.des)] = spd
            way_cnt += 1
            if way_cnt % 10_000 == 0
                @debug "processing $(way_cnt)-th way"
            end
        end
        
    end

    add_junction!(new_net)
    return new_net
end

"""
Create a new network by a induced graph and vmap from the function call of induced_subgraph
"""
function induced_network(net0::Network, sg, vmap)


    nodesdata = [net0.nodesdata[vmap[i]] for i in 1:nv(sg)]
    ## Change the index of the nodesdata 
    for (i, nd) in enumerate(nodesdata)
        nw_nd = update_idx(nd, i)
        nodesdata[i] = nw_nd
    end


    waydata = Dict(
       (e.src, e.dst) =>  net0.waydata[(vmap[e.src], vmap[e.dst])] for e in edges(sg)
    )

    for (key,way) in waydata
        new_way = update_src_des(way, key[1], key[2])
        waydata[key] = new_way
    end

    speeddata = Dict(
       (e.src, e.dst) => net0.speeddata[(vmap[e.src], vmap[e.dst])] for e in edges(sg)
    )

    ## Need to be careful about this...
    
    id2idx = Dict(
        nd.id => i for (i, nd) in enumerate(nodesdata)    
        # net0.id2idx[vmap[i]] => i for i in 1:nv(sg)
    )


    net = Network(;
        graph = sg,
        nodesdata = nodesdata,
        waydata = waydata,
        speeddata = speeddata,
        id2idx = id2idx,
        region = net0.region
    )

    add_junction!(net)

    return net

end

"""
Create a new network by removing the nodes that are isolated
"""
function induced_network_by_node(net0::Network, min_node_in_cc::Int)
    # max_com = largest_connected_component(net0.graph)
    commonents = weakly_connected_components(net0.graph)
    vlist = Int[]
    for com in commonents
        if length(com) > min_node_in_cc
            vlist = vcat(vlist, com)
        end
    end

    sg, vmap = induced_subgraph(net0.graph, vlist)

    return induced_network(net0, sg, vmap)
    
end

"""
Create a new network by removing the nodes that are isolated
"""
function induced_network_by_edge(net0::Network)
    elist = collect(edges(net0.graph))

    sg, vmap = induced_subgraph(net0.graph, elist)

    return induced_network(net0, sg, vmap)
end

"""
Find the "neighbors" that within one hop in node_set
"""
function find_neigh_bfs(net::Network, s::Int, node_set)
    q = Queue{Int}()

    n_node = nv(net)
    enqueue!(q, s)
    vlist = Int[]
    visited = falses(n_node)
    visited[s] = true
    while !isempty(q)
        u = dequeue!(q)
        visited[u] = true
        if (u!=s) && (u in node_set) && !(u in vlist)
            push!(vlist, u)
            continue
        end

        for v in outneighbors(net, u)
            if !visited[v]
                visited[v] = true
                enqueue!(q, v)
            end
        end

    end

    return vlist
end


function merge_one_degree_node(net0::Network; min_dis::Real = 1000.0, grade_tol::Real = 2.5e-3, region::String) 
    @info "Merging one degree nodes with min_dis=$(min_dis) and grade_tol=$(grade_tol)"


    ## We first create a new network by removing the nodes that are isolated
    # net = induced_network_by_node(net0, 10)

    # net = copy(net0)
    # max_com = largest_connected_component(net0.graph, false)
    max_com = net0.max_com
    sg, vmap = induced_subgraph(net0.graph, max_com)
    net = induced_network(net0, sg, vmap)
    # merge all one degree nodes
    
    for v in 1:nv(net)
    # for v in sort(net0.max_com)
        if v % 100_000 == 0
            perc = round(v / nv(net) * 100, digits=2)
            @debug "processing $(v)-th node, total $(nv(net)), current edges $(ne(net)) progress $(perc)%"
        end

        ## If the node is isolated or an endpoint, we do not need to merge it
        if indegree(net, v) == 0 || outdegree(net, v) == 0
            continue
        end
        
        ## u -> v -> w
        u = inneighbors(net, v)[1]
        w = outneighbors(net, v)[1]
        dis_in = distance(net, u, v)
        grade_in = grade(net, u, v)
        dis_out = distance(net, v, w)
        grade_out = grade(net, v, w)
        if (indegree(net, v) == 1) && (outdegree(net, v) == 1)
            if dis_in < min_dis || dis_out < min_dis || is_adj_edge_same_grade_level(net, v, grade_tol) || 
                !isvalid_grade(grade_in, g_min_theta, g_max_theta) || 
                !isvalid_grade(grade_out, g_min_theta, g_max_theta)
                    merge_one_degree_node!(net, v)
            end
        elseif indegree(net, v) == 1
            if dis_in < min_dis || !isvalid_grade(grade_in, g_min_theta, g_max_theta)
                # @debug "merging one indegree node $(v) with $(u) and $(w)"
                # merge_one_indegree_node!(net, v)
            end
        elseif outdegree(net, v) == 1
            if dis_out < min_dis || !isvalid_grade(grade_out, g_min_theta, g_max_theta)
                # merge_one_outdegree_node!(net, v)
            end
        end
    end

    @debug "Merging done, now removing isolated nodes"
    net = induced_network_by_node(net, 100)
    # net = induced_network_by_edge(net)

    net.region = region

    return net
end

function isvalid_grade(g::Real, min_g::Real, max_g::Real)
    if isnan(g)
        return false
    end
    if g < min_g || g > max_g
        return false
    end

    return true
end


"""
Check if the in and out edges are the same grade level
    
u -> v -> w
"""
function is_adj_edge_same_grade_level(net::Network, v::Int, tol::Real) 
    w = outneighbors(net, v)[1]
    u = inneighbors(net, v)[1]

    g1 = grade(net, u, v)
    g2 = grade(net, v, w)
    flag = abs(g1 - g2) < tol

    return flag
end

"""
Check if the node is too close to its neighbors.

This is used to check if we want to merge the node with its neighbors.
"""
function is_adj_edge_too_short(net::Network, node_idx::Int, min_dis::Real) 
    for v in outneighbors(net, node_idx)
        dis = distance(net, node_idx, v)
        if dis < min_dis
            return true
        end
    end

    for v in inneighbors(net, node_idx)
        dis = distance(net, v, node_idx)
        if dis < min_dis
            return true
        end
    end
    
    return false
end 

# is_isolated_node(net::Network, node_idx::Int) = outdegree(net, node_idx) == 0 && indegree(net, node_idx) == 0

# is_one_degree_node(net::Network, node_idx::Int) = outdegree(net, node_idx) == 1 && indegree(net, node_idx) == 1
"""
Merge a one indegree node with its neighbors

u -> v -> w1, & [, v -> w2]

The remaining graph only contains u -> w, and v becomes a isolated node.

min_dis is the minimum distance between the two nodes
"""
function merge_one_indegree_node!(net::Network, v::Int) 
    if indegree(net, v) != 1
        @debug "$(v) is not a one indegree node" indegree(net, v) outdegree(net, v)
        return
    end

    u = inneighbors(net, v)[1]

    out_neighbors_vec = collect(outneighbors(net, v))
    for w in out_neighbors_vec
        w == u && continue
        speed1 = net.speeddata[(u, v)]
        speed2 = net.speeddata[(v, w)]
        speed = (speed1 + speed2) / 2.0
        remove_way!(net, v, w)
        if !has_way(net, u, w) && (u != w)
            # add the way
            # @debug "adding way $(u) -> $(w)"
            addway!(net, u, w)
            net.speeddata[(u, w)] = speed
        end
    end

    remove_way!(net, u, v)

    if length(out_neighbors_vec) > 1
        # @debug "merging u -> v -> [w1, w2, w3]" Graphs.degree(net.graph, v) outdegree(net, u) indegree(net, u) length(out_neighbors_vec)
    end
end

"""
Merge a one outdegree node with its neighbors

[u1, u2] -> v -> w

The remaining graph only contains u -> w, and v becomes a isolated node.

min_dis is the minimum distance between the two nodes
"""
function merge_one_outdegree_node!(net::Network, v::Int) 
    if outdegree(net, v) != 1
        @debug "$(v) is not a one indegree node" indegree(net, v) outdegree(net, v)
        return
    end

    w = outneighbors(net, v)[1]

    for u in inneighbors(net, v)
        w == u && continue
        speed1 = net.speeddata[(u, v)]
        speed2 = net.speeddata[(v, w)]
        speed = (speed1 + speed2) / 2.0
        remove_way!(net, u, v)
        if !has_way(net, u, w) && (u != w)
            # add the way
            # @debug "adding way $(u) -> $(w)"
            addway!(net, u, w)
            net.speeddata[(u, w)] = speed
        end
    end

    remove_way!(net, v, w)
end



"""
Merge a one degree node with its neighbors

u -> v -> w

The remaining graph only contains u -> w, and v becomes a isolated node.

min_dis is the minimum distance between the two nodes
"""
function merge_one_degree_node!(net::Network, v::Int) 
    if indegree(net, v) != 1 || outdegree(net, v) != 1
        @debug "$(v) is not a one degree node" indegree(net, v) outdegree(net, v)
        return
    end

    return merge_one_indegree_node!(net, v)
end

function remove_way!(net, u::Int, v::Int)
    if !has_way(net, u, v)
        @warn "no way between $(u) and $(v)"
        return
    end

    # remove the way
    delete!(net.waydata, (u, v))
    delete!(net.speeddata, (u, v))
    rem_edge!(net.graph, u, v)
end