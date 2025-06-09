"""
"""


"""
since we focus on largest component (net.max_com), 
we need to find the closest node in the junctions set to avoid
routing between two disconnected components
"""
function closest_node_in_max_com(net::Network, n::Int64)::Int64
    if n in net.max_com 
        return n
    end
    dists = [distance2d(net, n, n1) for n1 in net.max_com]
    (val,idx) = findmin(dists)
    return net.max_com[idx]
end

closest_node_in_max_com(net::Network, latlon) = closest_node_in_max_com(net, closest_node(net, latlon))

"""
Return the closest node in the network to the given latlon, return the idx
"""
function closest_node(net::Network, latlon)::Int64
    dists = (distance2d(latlon, LatLon(node)) for node in net.nodesdata)
    (val,idx) = findmin(dists)
    # @show val idx
    node_idx = net.nodesdata[idx].idx
    return node_idx
end

"""
find closest_node for a vector of nodes with KDTree
"""
function closest_node(net::Network, latlons::Vector{T}; max_com_flag::Bool = false) where T <: LatLon
    if max_com_flag
        data = [[getnode(net, i).lat, getnode(net, i).lon] for i in net.max_com]
    else
        data = [[n.lat, n.lon] for n in net.nodesdata]
    end
    data = reduce(hcat, data)
    tree = NN.KDTree(data)
    node_idx_vec = [
        NN.nn(tree, [x.lat, x.lon])[1] for x in latlons
    ]
    if max_com_flag
        node_idx_vec = [net.max_com[i] for i in node_idx_vec]
    end
    return node_idx_vec
end
