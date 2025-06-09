"""
The nework part that store the graph representation of the highway system
"""

#const k_region = "alabama"
const g_net_data_path = normpath(joinpath(k_root_dir,"data","map","osm-data"))
global g_bin_net_data_dir = normpath(joinpath(k_root_dir,"data","map","bin"))

"""
The Graph/MetaGraph lib does not provide the option to use dictionary.
This is a wrapper for it.

graph: the graph that store the basic information

id2idx: The dictionary that map from id to index

junctions: A junction in the setting of cfo is a node with more than two neighbors.

max_com: the maximum connected component

waydata: a dictionary that map a pair of int to an way obj

region: 

speeddata:

meta: store the information of other meta data
"""
@with_kw mutable struct Network{NT <: AbstractNode, WT <: AbstractWay} <: AbsNet
	graph::DiGraph{Int64} = DiGraph()
	id2idx::Dict{Int64,Int64} = Dict()  # store the information from 
	junctions::Vector{Int64} = []
	max_com::Vector{Int64} = [] # maximum connected component
	nodesdata::Vector{NT} = StaticNode[] # store the nodes with idx as a key
	waydata::Dict{Tuple{Int64, Int64}, WT} = Dict{Tuple{Int, Int}, StaticEdge}()
	region::String
	# speeddata::OrderedDict{Tuple{Int64, Int64}, Float64} = OrderedDict()
	speeddata::Dict{Tuple{Int64, Int64}, Float64} = Dict()
    meta::Dict{Symbol, Any} = Dict() # store the information of other meta data
end
Network(region::String) = Network(;region=region)
Graphs.edges(net::Network) = Graphs.edges(net.graph)
Graphs.inneighbors(net::Network, nd::Int64) = inneighbors(net.graph, nd)
Graphs.outneighbors(net::Network, nd::Int64) = outneighbors(net.graph, nd)
Graphs.indegree(net::Network, nd::Int64) = indegree(net.graph, nd)
Graphs.outdegree(net::Network, nd::Int64) = outdegree(net.graph, nd)
Graphs.all_neighbors(net::Network, nd::Int64) = all_neighbors(net.graph, nd)

function Base.getproperty(net::Network, name::Symbol)
    if name in fieldnames(Network)
        return Base.getfield(net, name)
    elseif name == :g
        return Base.getfield(net, :graph)
    elseif haskey(net.meta, name)
        return net.meta[name]
    else
        throw(ErrorException("Network has no property $name"))
    end
end

function Base.setproperty!(net::Network, name::Symbol, x)
    if name in fieldnames(Network)
        return Base.setfield!(net, name, x)
    else
        return net.meta[name] = x
    end
end

function Base.copy(net0::Network)
    net = Network(;
        region = net0.region,
        nodesdata = net0.nodesdata,
        waydata = copy(net0.waydata),
        speeddata = copy(net0.speeddata),
        junctions = copy(net0.junctions),
        max_com = copy(net0.max_com),
        graph = copy(net0.graph),
        id2idx = copy(net0.id2idx),
        meta = copy(net0.meta),
    )
    return net
end


function check_self_loop(net::Network)
    for e in edges(net.graph) 
        if e.src == e.dst
            @warn "checking self loop, src==des==$(e.src) !"
        end
    end
end

function addnode!(net::Network{NT}, node::NT) where NT <: AbstractNode
	add_vertex!(net.graph)
	idx = node.idx
	@assert(idx == Graphs.nv(net.graph))
    @assert(idx == length(net.nodesdata)+1)
	push!(net.nodesdata, node)
	net.id2idx[node.id] = node.idx
end

function has_node(net::Network, idx::Int)
    idx <= length(net.nodesdata) 
end

function has_way(net::Network, u::Int, v::Int)
    return haskey(net.waydata, (u,v)) 
end

"""
"""
function get_tag_dict(node::EzXML.Node)::Dict
	#@show(typeof(node))
	tag_dict = Dict()
	for tag in eachelement(node)
		if nodename(tag) == "tag"
			k,v = tag["k"], tag["v"]
			tag_dict[k] = v
		end
	end
	return tag_dict
end

region2osmfile(region::String) = joinpath(g_net_data_path,region,"$(region)-highway.osm")
region2nodefile(region::String) = joinpath(g_net_data_path,region,"$(region)-highway-node-ele.json")
region2wayfile(region::String)::String = joinpath(g_net_data_path,region,"$(region)-highway-way.json")
region2nodecsv(region::String) = joinpath(g_net_data_path, region, "$(region)_node.csv")
region2edgecsv(region::String) = joinpath(g_net_data_path, region, "$(region)_edge.csv")

id2idx(net::Network, id::Integer) = net.id2idx[id]
id2idx(net::Network, id::String) = id2idx(net, parse(Int64, id))
    
# function id2idx(net::Network, id)
# 	if isa(id, String)
# 		id = parse(Int64, id)
# 	end
# 	return net.id2idx[id]
# end
idx2id(net::AbsNet, idx::Int) = getnode(net, idx).id

function printnetstat(net::Network)
	n_nodes = nv(net)
	nodes_mem = Base.summarysize(net.nodesdata)/1024/1024
	mem_per_node = nodes_mem / n_nodes

	n_ways = nway(net)
	edge_mem = Base.summarysize(net.waydata)/1024/1024
	mem_per_edge = edge_mem / n_ways
	@info "graph used $(Base.summarysize(net.graph)/1e6) M"
	@info "nodesdata used $(nodes_mem) M, $(mem_per_node*1024*1024) Byte per node"
	@info "waydata used $(edge_mem) M, $(mem_per_edge*1024*1024) Byte per edge"
	n_junections = n_junc(net)
	n_edges = ne(net)
	@info "$(n_junections) junctions in net."
	@info "$(n_nodes) nodes in net"
	@info "$(n_ways) ways in net"
	@info "$n_edges edges in net"
end

function distance(net::Network, idx1::Int64, idx2::Int64)::Float64
    # return distance3d(net, idx1, idx2)
    way = getway(net, idx1, idx2)
    # @assert(length(way.nodes) == 2)
    val::Float64 = way.distances[1]
    return val
end

function grade(net::Network, idx1::Int64, idx2::Int64)
    way = getway(net, idx1, idx2)
    # @assert(length(way.nodes) == 2)
    return way.grades[1]
    #g1 = net.grade_data[idx1, idx2]
    #if g1 == 0
    #    g = grade(net.nodesdata[idx1], net.nodesdata[idx2])
    #    net.grade_data[idx1, idx2] = g
    #    return g
    #else
    #    return g1
    #end
end

distance2d(net::Network, idx1::Int64, idx2::Int64) = distance2d(net.nodesdata[idx1], net.nodesdata[idx2])
distance3d(net::Network, idx1::Int64, idx2::Int64) = distance3d(net.nodesdata[idx1], net.nodesdata[idx2])


getway(net::Network, idx1::Int, idx2::Int)= net.waydata[(idx1, idx2)]	
getway(net::Network, e::Graphs.Edge) = getway(net, e.src, e.dst)
islink(net::Network, e::Graphs.Edge) = islink(getway(net, e))



getnode(net::Network, idx::Int64) = net.nodesdata[idx]	
function getnode(net::Network{NT, WT}, idx_v::Vector{Int64}) where {NT <: AbstractNode, WT <: AbstractWay}
    n = length(idx_v)
    node_v = Vector{NT}(undef, n)
    for i = 1:n
        node_v[i] = getnode(net, idx_v[i])
    end
    return node_v
end
isjunction(net::Network, idx::Int64) = (idx in net.junctions)
n_junc(net::Network) = length(net.junctions)


nway(net::Network) = length(net.waydata)
Graphs.nv(net::Network) = length(net.nodesdata)
Graphs.ne(net::Network) = sum(length(e.nodes)-1 for e in values(net.waydata))

function get_max_speed(net::Network, i::Int64, j::Int64)::Float64
    data_dict = net.speeddata
    @inbounds data::Float64 = data_dict[(i, j)]
    return data == 0.0 ? g_max_highway_speed : data
end
get_max_speed(net::Network, way::Way) = get_max_speed(net, way.src, way.des)




