"""
"""

"""
Single source shortest to all other vertices.
return the distance vec from src to each node 
"""
function shortest_path_all(net::AbsNet, src::Int64, distmx)

    n_node::Int = nv(net)

	# open_set = PriorityQueue{Int, Float64}()
	open_set = HeapPriorityQueue{Int, Float64}()
    enqueue!(open_set, src, 0.0)

    # Delete the visited vector to allow enqueue again: useful when there is negative weight
	count = zeros(Int64, n_node)

    g_score = fill(Inf, n_node)
    g_score[src] = 0

    inqueue = falses(n_node,1)
    inqueue[src] = true

	parents = zeros(Int64, n_node)

    #visited(u::Int) = visitonce && count[u] == 0

	cnt = 0
	while !isempty(open_set)
        # @debug open_set
		u::Int = dequeue!(open_set)
        count[u] += 1
        inqueue[u] = false
		# @debug "visting $u"
        neigh_vec = outneighbors(net, u) 
        iv::Int = 1
        # @inbounds for iv in 1:length(neigh_vec)
        ## Note: for loop can be slow and critical, so we manually write a while loop here.
        # while iv <= length(neigh_vec)
        #     v::Int = neigh_vec[iv]
		@inbounds for v::Int in neigh_vec
            # This is tricky, to avoid negative loop
            # if visitonce && count[v] != 0
            #     #@debug "node $v visited with count $(count[v]), skip it."
            #     continue
            # end

            # dist = getdist(net, u, v)[1]
            @inbounds dist::Float64 = distmx[u,v]
            cnt += 1
            tentative_g_score = g_score[u] + dist
            if tentative_g_score < g_score[v]
                g_score[v] = tentative_g_score
                parents[v] = u
                priority = tentative_g_score 
                # priority = tentative_g_score + heuristic(v)
                enqueue!(open_set, v, priority)
                if !inqueue[v]
                    inqueue[v] = true
                end
            end
            if count[v] >= n_node
                # print_sparse(distmx, :distmx)
                # @debug g_score
                @warn "There is a negative Cycle in the network. Total $cnt getdist called." maxlog=10
                throw(Graphs.NegativeCycleError())     
                # return Vector{Int}()
                return g_score
            end
            iv += 1
		end ## end of iteration for neigh_vec
	end

	# @debug "Total $cnt getdist called."

    return g_score

end