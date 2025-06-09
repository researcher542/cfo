"""
Get the augmented weight for each edge in the graph.
"""

"""
Get the distance matrix for each pair of charging stations.
This is the one that should be called by others.
mx[i][u,v] is the SP^i(u,v)
"""
function ds_get_cs_sp_mx(prob::AbsCfoProb, cs_vec::Vector{ChargeStation}, dual::DsDual, N::Int, option::DsOption{DsNetwork}) 
    (;net, src, des, ev) = prob
    # The low-level distance matrix 
    @timeit g_to "ds_get_dist_mx" low_dist_mx = ds_get_dist_mx(net, dual, ev, N, option, prob.objtype)
    # append src,des in the node list
    cs_idx_vec = vcat(src, des, [cs.idx for cs in cs_vec])
    unique!(cs_idx_vec)
    I = [(i, j)[1] for i in cs_idx_vec for j in cs_idx_vec]
    J = [(i, j)[2] for i in cs_idx_vec for j in cs_idx_vec]
    V = zeros(length(I)) 
    n_node = nv(net)
    cs_dist_mx = [spzeros(n_node, n_node) for _ in 1:N+1]
    for i in 1:N+1
        # @debug i
        fill!(V, Inf)
        # getdist(net, u, v) = low_dist_mx[i][u,v]
        # iv::Int = 1
        GC.gc(false)
        Threads.@threads :greedy for isrc in eachindex(cs_idx_vec)
            # @timeit g_to "shortest path" dist_mx_i = shortest_path(Astar(), net, src, src; breakearly=false, outputdists=true, getdist=getdist)
            src1 = cs_idx_vec[isrc]
            ## The most time-consuming and most memory-intensive part.
            # @timeit g_to "shortest_path_all" 
            dist_mx_i = shortest_path_all(net, src1, low_dist_mx[i])
            # @timeit g_to "shortest path" ds = Graphs.dijkstra_shortest_paths(net.g, src, low_dist_mx[i])
            # dist_mx_i = ds.dists
            # @timeit g_to "shortest path" dist_mx_i = Graphs.spfa_shortest_paths(net.g, src, low_dist_mx[i])
            for (ides,des1) in enumerate(cs_idx_vec)
                ivv = (isrc-1) * (length(cs_idx_vec)) + ides
                V[ivv] = dist_mx_i[des1]
                # cs_dist_mx[i][src, des] = dist_mx_i[des]
                # @show isrc ides iv ivv
                # @assert iv == ivv 
                # iv += 1
            end
        end
        AA = sparse(I, J, V, n_node, n_node)
        cs_dist_mx[i] = AA
        # @assert cs_dist_mx[i] == AA
    end
    return cs_dist_mx
end



"""
Get the distance matrix on the transportation graph
Note that this function is a aux function for computing augmented cost between two charging stations.

return:  dist_mx[istage][u, v] = cost
"""
function ds_get_dist_mx(net::Network, dual::DsDual, ev::EV, N::Int, option::DsOption, objtype::AbstractObjective)
    n_node = nv(net)
    n_edge = ne(net)
    dist_mx = [spzeros(n_node, n_node) for _ in 1:N+1]
    I = [e.src for e in edges(net.g)]
    J = [e.dst for e in edges(net.g)]
    V = zeros(n_edge)
    # dist_mx = [RobinDict{Tuple{Int, Int}} for _ in 1:N+1]
    edge_vec = collect(edges(net.g))
    for i in 1:N+1
        V = zeros(n_edge)
        Threads.@threads :greedy for ie::Int in 1:length(edge_vec)
            e = edge_vec[ie]
            u,v = e.src, e.dst
            # @timeit g_to "ds_get_dist" 
            dist = ds_get_dist(net, i, u, v, dual, ev, option, objtype)[1]
            # dist_mx[i][u,v] = dist
            V[ie] = dist
            # ie += 1
        end
        AA = sparse(I, J, V)
        # @assert dist_mx[i] == AA
        dist_mx[i] = AA
    end
    return dist_mx
end

"""
Get the augmented distance for subpath i between node u and v,
    return: (cost, t)
"""
function ds_get_dist(net::Network, i::Int, u::Int, v::Int, dual::DsDual, ev::EV, option::DsOption, objtype::AbstractObjective)::Tuple{Float64, Float64}
    μi = dual.μ[i]
    λi = dual.λ[i]
    min_t, max_t = get_minmax_t(net, u, v, ev)
    (;t_mul,b_mul) = option
    if objtype == ObjTime()
        λi += 1.0
    end
    if λi == 0 && μi == 0
        # t::Float64 = (min_t+max_t)/2
        t::Float64 = min_t
        cost::Float64 = 0.0
    elseif μi == 0
        @assert(λi > 0)
        t = min_t
        cost = λi*t*t_mul
    else
        # if objtype == ObjTime()
        #     λ0 = ( (λi + 1) *t_mul) / (μi*b_mul)
        #     t, opt_spd = augmentedtime(net, u, v, λ0, ev)
        #     energy_cost = energy_cost_on_road(net, u, v, ev, t)
        #     cost = energy_cost*μi*b_mul + (λi+1)*t*t_mul
        # else
        # end
        λ0 = (λi*t_mul) / (μi*b_mul)
        t, opt_spd = augmentedtime(net, u, v, λ0, ev)
        energy_cost = energy_cost_on_road(net, u, v, ev, t)
        cost = energy_cost*μi*b_mul + λi*t*t_mul

        ## TODO: This is tricky... to avoid negative weight
        if option.pos_cost && energy_cost < 0.0 && cost < 0
            cost = 0.0
        end
    end

    return (cost, t)
end

