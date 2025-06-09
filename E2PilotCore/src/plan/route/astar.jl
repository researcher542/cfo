
using FastPriorityQueues

function get_heuristic(net::Network, src::Int, des::Int, veh::AbstractVehicle, objtype::AbstractObjective; visitonce::Bool=false)
    function getdist(net, u, v, lambda)
        dis = distance(net, u, v)
        theta = grade(net, u, v)
        max_spd = get_max_speed(net, u, v, veh)
        # time = distance / (spd+1e-4)
        opt_spd = optimal_speed(veh, theta, 0.0, g_min_highway_speed,max_spd)
        e = output_power(veh, theta, opt_spd) * (dis / opt_spd)
        # @info  "min_time called $time"
        return e
    end
    dists = shortest_path(Astar(), net, src, des; breakearly=false, outputdists=true, getdist=getdist, visitonce=visitonce)
    hd = dists[des]
    return x-> (hd- dists[x])
end


"""
Shortest path algorithm with A^Star algorithm. 
Note that if heurstic zero, then it is simply a Dijkstra.
Also, if breakearly is false, then it is the SPFA algorithm (?need to check)
Refer to https://github.com/JuliaGraphs/Graphs.jl/blob/master/src/shortestpaths/astar.jl
and http://en.wikipedia.org/wiki/A%2A_search_algorithm

getdist: A function (network, u, v, lambda)

heurstic: A heurstic used in A star, with g(u) be the heurstic score for node u. Note that if satisfies g(u)-g(v) < dist(u,v), then the result is still optimal.

return:

    path: Vector{Int}
"""
function shortest_path(type::Astar, net::AbsNet, src::Int64, des::Int64 ; 
        getdist::Function=(net, i, j, lambda)->distance(net, i, j), 
        lambda = 0.0,  heuristic::Function=x->0.0, 
        breakearly::Bool=false, 
        outputdists::Bool=false, 
        visitonce::Bool=false, 
        distmx = nothing, getdistargs = (),
        throw_err_flag::Bool = true,
        )

    if !isnothing(distmx)
        getdist = (net, u, v, lam) -> distmx[u, v]
    end
    n_node::Int = nv(net)

	# open_set = PriorityQueue{Int, Float64}()
	open_set = HeapPriorityQueue{Int64, Float64}()
    enqueue!(open_set, src, 0.0)

    # I = [e.src for e in edges(net.g)]
    # J = [e.dst for e in edges(net.g)]
    # V = fill(Inf, length(I))
    # distmx = sparse(I, J, V)
    # distmx = spzeros(n_node, n_node)

    # Delete the visited vector to allow enqueue again: useful when there is negative weight
    max_int = 2^31 - 1
    @assert (n_node <= max_int) && (n_node > 0)
	count = zeros(Int32, n_node)

    g_score = fill(Inf, n_node)
    g_score[src] = 0

    f_score = fill(Inf, n_node)
    f_score[src] = heuristic(src)

    inqueue = falses(n_node, 1)
    inqueue[src] = true

	parents = zeros(Int32, n_node)

    #visited(u::Int) = visitonce && count[u] == 0

	cnt = 0
	while !isempty(open_set)
        # @debug open_set
		u = dequeue!(open_set)
        if visitonce && count[u] != 0
            # @debug "node $u visited with count $(count[u]), skip it."
            continue
        end
        count[u] += 1
        inqueue[u] = false
		# @debug "visting $u"
		if breakearly && (u == des)
			break
		end
		for v in outneighbors(net, u) 
            # This is tricky, to avoid negative loop
            if visitonce && count[v] != 0
                #@debug "node $v visited with count $(count[v]), skip it."
                continue
            end

            if !isempty(getdistargs)
                dist = getdist(net, u, v, lambda, getdistargs...)[1]
            else
                dist = getdist(net, u, v, lambda)[1]
            end

            cnt += 1
            # @debug "For edge ($u, $v), weight is $(dist)"
            tentative_g_score = g_score[u] + dist
            if tentative_g_score < g_score[v]
                g_score[v] = tentative_g_score
                parents[v] = u
                priority = tentative_g_score + heuristic(v)
                # open_set[v] = priority
                enqueue!(open_set, v, priority)
                if !inqueue[v]
                    inqueue[v] = true
                end
            end
            if count[v] >= n_node
                # print_sparse(distmx, :distmx)
                # @debug g_score
                @warn "There is a negative Cycle in the network"
	            @debug "Total $cnt getdist called."
                # @exfiltrate
                throw(Graphs.NegativeCycleError())     
                return Vector{Int}()
            end
            
		end
	end
	# @debug "Total $cnt getdist called."
	if parents[des] == 0 && !(outputdists)
		@warn "Cannot find a path from $src to $des"
        if throw_err_flag
            @exfiltrate
            throw(ErrorException("Cannot find a path from $src to $des"))
        end
        return Vector{Int}()
	end

    if outputdists
        return g_score
    end

    path = retrieve_path(src, des, parents)
    return path
end

function retrieve_path(src::Int, des::Int, parents::Vector) 
    path::Vector{Int64} = [des]
	p = parents[des]
	while p != src
        if p in path
            @error "There is a loop in the path! p=$p already in path=$path"
            @assert false
        end
		pushfirst!(path, p)
		p = parents[p]
	end
	pushfirst!(path, p)
    return path 
end