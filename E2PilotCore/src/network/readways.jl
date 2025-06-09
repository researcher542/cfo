"""
"""

"""
Note: The osm data format, ways can be something like the following. And osm treated it as two ways but not three ways....

_______________
____    _______
	|  |
	|  |
	|  |
To tackle this, we first add edges within the way and merge them later.
"""
function readways_osm!(net::Network, region::String = "test")
    @info "reading ways from osm file $(region)"
	function addjunction(net::Network, tag_dict::Dict, node::Int64)
		k = "highway"
		if k in keys(tag_dict)
			v = tag_dict[k]
			if occursin("link", v)
				push!(net.junctions, node)	
				return true
			end
		end
		return false
	end
	open(EzXML.StreamReader, region2osmfile(region)) do reader
		idx = 0
		for typ in reader
			if reader.name == "way" && typ == EzXML.READER_ELEMENT
				## get node list
				way = expandtree(reader)
				node_vec = []
				tag_dict = Dict()
				idx += 1
				if idx % 10000 == 0
					@debug "processing $(idx)-th ways"
				end
				for xml_node in eachelement(way)
					if nodename(xml_node) == "nd"
						id = str2int(xml_node["ref"])
                        if !haskey(net.id2idx, id)
                            @warn "id $id not in net.id2idx"
                            continue
                        end
						push!(node_vec, id2idx(net, id))
					elseif nodename(xml_node) == "tag"
						k,v = xml_node["k"], xml_node["v"]
						tag_dict[k] = v
					end
				end

				src = node_vec[1]
				des = node_vec[end]
				for i in 1:length(node_vec)-1
					src = node_vec[i]
					des = node_vec[i+1]
					#edge = Way(src, des, str2int(way["id"]), [src, des])
					#net.waydata[(src,des)] = edge
					add_edge!(net.graph, src, des)
				end
				#edge = Way(src, des, str2int(way["id"]), false, node_vec)
				islink = addjunction(net, tag_dict, src)
				#edge.islink = islink
				#net.waydata[(src,des)] = edge
				#@debug values(tag_dict)
				#net.waydata[(src,des)] = edge
			end
		end
	end
	return net
end

function largest_connected_component(g::AbstractGraph, strongly::Bool = true)
    if strongly
	    components = strongly_connected_components(g)
    else
	    components = weakly_connected_components(g)
    end
	(_, idx) = findmax(length, components)
	max_com = components[idx]
	return max_com
end
largest_connected_component(net::Network) = largest_connected_component(net.graph)


"""
Merge the edges that only has 1-out-degree and 1-out-degree in the network, and create the ways.
Only do this for the largest connected components.
    !!! Cannot deal with loop, depriated for now
"""
function merge!(net::Network)
	"""
	Add a way given source node
	"""
	function addway_src!(net::Network, src::Int64)
		#@debug "calling addway_src!"
		des_vec = []
        neigh = outneighbors(net, src)
		for nd in neigh
			node_vec = [src, nd]
			out_nodes = outneighbors(net, nd)
			nd1 = nd
			while !isfork(net, nd1) && !is_endpoint(net, nd1)
				nd1 = out_nodes[1]
				push!(node_vec, nd1)
				out_nodes = outneighbors(net, nd1)
			end
			des = nd1
			push!(des_vec, des)
			islink_ = isjunction(net, src)
            if src == des
                @warn "src==des==$src !"
                @assert false
                #throw(InitError(:test))
            end
			way = Way(net, src, des, -1, islink_, node_vec)
			addway!(net, way)
		end
	end

	#max_com = largest_connected_component(net)	
	for nd in vertices(net.graph)
		if isfork(net, nd)
			addway_src!(net, nd)
		end
	end

	# clean the minor edges
	for e in edges(net)
		if !haskey(net.waydata, (e.src, e.dst))
			rem_edge!(net.graph, e)
		end
	end

end

isfork(net::Network, nd::Int64)  = (outdegree(net,nd) > 1) || (indegree(net, nd) > 1)
function is_endpoint(net::Network, nd::Int64)
    return length(all_neighbors(net, nd)) <= 1
end

"""
simply add all ways to the network
"""
function addways!(net::Network)
    for e in edges(net.graph)
        src = e.src
        des = e.dst
	    # way = Way(net, src, des, -1, false, [src, des])
        addway!(net, src, des)
    end
end

function addway!(net::Network{NT, NW}, src::Int, des::Int) where {NT, NW}
    way = NW(net, src, des, -1, false, SA[src, des]) 
    addway!(net, way)
    if !haskey(net.speeddata, (src, des))
        net.speeddata[(src, des)] = 0.0
    end
end

"""
"""
function addway!(net::Network, way::StaticEdge) 
	net.waydata[(way.src, way.des)] = way
	if !has_edge(net.graph, way.src, way.des)
		add_edge!(net.graph, way.src, way.des)
	end
end


"""
"""
function addway!(net::Network, way::Way) 
	net.waydata[(way.src, way.des)] = way
	nodes = way.nodes
	# remove the intermediate edges
	if length(nodes) >= 3
		for i in 1:length(nodes)-1
			rem_edge!(net.graph, nodes[i], nodes[i+1])
		end
	end
	if !has_edge(net.graph, way.src, way.des)
		add_edge!(net.graph, way.src, way.des)
	end
end


reduce!(net::Network, e::StaticEdge) = nothing

"""
reduce the nodes in the edge by its distance and its grade.
"""
function reduce!(net::Network, e::AbstractWay)
	nodes = getnode.((net,), e.nodes)
	nodes1 = reduce_by_distance(nodes)
	nodes2 = reduce_by_grade(nodes1)
	e.nodes = [n.idx for n in nodes2]
	return e
end

function reduce_edge!(net::Network)
    n_edge1 = ne(net) 
	for edge in values(net.waydata)	
		reduce!(net, edge)
	end
    n_edge2 = ne(net) 
    @debug "reduce edge from $n_edge1 to $n_edge2"
end