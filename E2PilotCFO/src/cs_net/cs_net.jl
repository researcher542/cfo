

"""
t_cost_vec: store a vector of (time, cost) pairs
"""
@with_kw struct CsEdge{TT <: Real}
    min_cost::TT = Inf
    max_cost::TT = Inf
    min_t::TT = Inf
    max_t::TT = Inf
    src::Int = -1
    des::Int = -1
    dis::TT = Inf
    # t_cost_vec::Vector{Tuple{TT,TT}} = Tuple{Float64,Float64}[]
    t_vec::Vector{TT} = Float64[]
    c_vec::Vector{TT} = Float64[]
end

function (edge::CsEdge{TT})(t::TT, ext::Bool = false)::TT where TT <: Real
    (;t_vec, c_vec) = edge
    # t_vec = [x[1] for x in edge.t_cost_vec] 
    # cost_vec = [x[2] for x in edge.t_cost_vec] 
    if length(t_vec) == 1
        return c_vec[1]
    else
        @assert issorted(t_vec) 
        return ep.interpolate_pwl(t, t_vec, c_vec, ext)
        # pwl = PWL(; x=t_vec, y=c_vec)
        # return pwl(t, ext)
    end
end

"""
This can be tricky.
Here, id means the idx from original graph, not their id...
"""
@with_kw mutable struct CsNet{T,G} <: AbsNet
    # data matrix that stores the edge data
    data::T
    g::G = DiGraph()
    id2idx::Dict{Int,Int} = Dict()
    idx2id::Dict{Int,Int} = Dict()
    cs_nei_flag::Bool = false # if true, then we need to add an edge only if there is a path not going through other charging stations.
end

Base.isempty(cs_net::CsNet) = Base.isempty(cs_net.data)

Graphs.nv(net::CsNet) = maximum(keys(net.data))
Graphs.ne(net::CsNet) = sum(map(length, values(net.data)))
function outneighbors(net::CsNet, u::Int)
    # If there is no key, then it is probably because the src or des is not a charging station.
    edges = net.data[u]
    return [e[1] for e in edges]
end

function inneighbors(cs_net::CsNet, v::Int)
    v_idx = cs_net.id2idx[v]
    in_nei = inneighbors(cs_net.g, v_idx)
    return [cs_net.idx2id[u_idx] for u_idx in in_nei]
end

function Graphs.has_edge(cs_net::CsNet, i::Int, j::Int)
    return has_edge_net(cs_net, i, j)
    # ii = cs_net.id2idx[i]
    # jj = cs_net.id2idx[j]
    # return Graphs.has_edge(cs_net.g, ii, jj)
end

function has_edge_net(cs_net::CsNet, i::Int, j::Int)
    if !haskey(cs_net.data, i)
        return false
    else
        edge = get_edge(cs_net, i, j)
        return !isnothing(edge)
        # des_vec = [e[1] for e in cs_net.data[i]]
        # if j in des_vec
        #     return true
        # end
    end
    return false
end

function get_edge(net::CsNet, i::Int, j::Int)
    edges = net.data[i]
    for e in edges
        if e[1] == j
            ee = e[2]
            if isempty(ee.t_vec)
                return nothing
            else
                return ee
            end
        end
    end
    return nothing
end

"""
"""
function Base.getindex(net::CsNet, i::Int, j::Int)
    return get_edge(net, i, j)
    # throw(ErrorException("edge not exists"))
end

# """
# Get one cs_net edge on the original graph.
# Should be careful about the index
# """
# function get_cs_net_edge(net::Network, i::Int, j::Int, ev::EV)
#     dis = distance(net, i, j)
#     return cs_net_get_one_edge(net, i, j, ev, dis) 
# end

function distance(net::CsNet, i::Int, j::Int) 
    e::CsEdge = get_edge(net, i, j)
    return e.dis
end

function CsNet(n::Int; kwargs...)
    data = Dict{Int, Vector{Tuple{Int64, CsEdge{Float64}}}}()
    return CsNet(; data=data, kwargs...)
end

isinf(x) = Base.isinf(x)
function isinf(e::CsEdge)
    for sym in fieldnames(CsEdge)
        val = getfield(e, sym)
        if isinf(val)
            return true
        end
    end
    return false
end

function add_edge_net!(net::CsNet, i::Int, j::Int, edge::CsEdge)
    data = (j, deepcopy(edge)) 
    if !haskey(net.data, i)
        net.data[i] = CsEdge[]
    end
    push!(net.data[i], data)
    return net
end

"""
fastflag: if true, we restrict the speed to maximum speed
"""
function minmax_t(net::CsNet, i::Int, j::Int, veh::AbstractVehicle, fastflag::Bool = false)
    if (i==j)
        return (0.0, 0.0)
    end
    edge = get_edge(net, i, j)
    if isnothing(edge)
        error("cs_net edge ($i, $j) not found.")
    end
    (;min_t, max_t) = edge
    if fastflag
        return (min_t, min_t)
    end
    # dis = distance(net, i, j)
    # min_spd,max_spd = get_minmax_speed(net, i, j, veh)
    # min_t, max_t = dis/max_spd, dis/min_spd
    return min_t, max_t
end

function closest_cs_node(net::Network, i::Int)
    min_dis = Inf
    idx = -1
    for cs in net.cs_vec
        dis = distance3d(net, cs.idx, i)
        if dis < min_dis
            min_dis = dis
            idx = cs.idx
        end
    end
    return idx
end

include("get_cs_net.jl")
include("path.jl")
include("save.jl")
include("var.jl")