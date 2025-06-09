

"""
reduce a network to a smaller one that is within a box given by latlons
"""
function reduce_net(net::Network{NT, WT}, region::String = "reduced_net"; min_lat::T = -Inf, max_lat::T = Inf, min_lon::T = -Inf, max_lon::T = Inf) where {T <: Float64, NT <: AbstractNode, WT <: AbstractWay}


    reduced_vlist = Int[]
    for node in net.nodesdata
        (;lat,lon) = node
        if (min_lat <= lat <= max_lat)  && (min_lon <= lon <= max_lon)
            push!(reduced_vlist, node.idx)
        end
    end

    sg, vmap = induced_subgraph(net.graph, reduced_vlist)

    net1 = induced_network(net, sg, vmap)
    net1.region = region
    return net1


    # function node_old2new(net, new_net, old_idx::Int)
    #     node = getnode(net, old_idx)
    #     id = node.id
    #     new_idx = get(new_net.id2idx, id, -1)
    #     return new_idx
    # end

    # new_net::Network = Network(region)
    # node_cnt = 1
   

    # # new_net.speeddata = deepcopy(net.speeddata)
    # for way in values(net.waydata)
    #     src = node_old2new(net, new_net, way.src)
    #     des = node_old2new(net, new_net, way.des)
    #     if src == -1 || des == -1
    #         continue
    #     end
    #     src_node_new = getnode(new_net, src)
    #     des_node_new = getnode(new_net, des)
    #     theta = grade(src_node_new, des_node_new)
    #     dis = distance3d(src_node_new, des_node_new)
    #     if has_node(new_net, src) && has_node(new_net, des)
    #         # new_way = Way(; src=src, des=des, id=way.id, nodes=[src, des], grades=[theta], distances=[dis])
    #         new_way = WT(new_net, src, des=des, id=way.id, nodes=[src, des], grades=[theta], distances=[dis])
    #         @assert new_way.grades == way.grades "$new_way $way"
    #         new_way.src = src; new_way.des = des;
    #         addway!(new_net, new_way)
    #         new_net.speeddata[(new_way.src, new_way.des)] = net.speeddata[(way.src, way.des)]
    #     end
    # end

    # add_junction!(new_net)
    return new_net 
end

