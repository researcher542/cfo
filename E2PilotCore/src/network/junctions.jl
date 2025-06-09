

"""
We define a junction is a node with more than two neighbors (>2) including in-neighbors and out-neighbors.
Add junction vector to the network structure.
Also add max_com vector to the network structure.
cs_net_flag: If we do this for the cs_net.
"""
function add_junction!(net::Network; min_nei::Int = 3, cs_net_flag::Bool = false)
	origin_n_junc = n_junc(net)
	#@debug "cleaning junctions with njunction=$(n_junc(net))"
    n = nv(net)
	max_com = largest_connected_component(net)	
	net.max_com = max_com
    junctions = Int[]
    if cs_net_flag
        cs_idx_r = [cs.idx_r for cs in net.cs_vec]
    end
    for inode in 1:n
        nei_vec = all_neighbors(net.graph, inode)
        n_nei = (length(nei_vec))
        if cs_net_flag && (n_nei <= min_nei+1)
            if (inode in cs_idx_r)
                continue
            end
            if !isempty(intersect(nei_vec, cs_idx_r))
                continue
            end
        end
        if n_nei >= min_nei
            push!(junctions, inode)
        end
    end
    net.junctions = junctions

	
    unique!(net.junctions)
	
	return net
end