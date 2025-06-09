

abstract type AbstractWay end

"""
A immutable struct for the way data.

We might need to avoid mutable structs inside

We assume there is no intermediate nodes in an edge
"""
@with_kw struct StaticEdge <: AbstractWay
    src::Int64
    des::Int64
    id::Int64
    grade::Float64
    distance::Float64
end

function StaticEdge(net, src::Int, des::Int, id, islink, nodes)
    n1 = getnode(net, src)
    n2 = getnode(net, des)
    g = grade(n1, n2)
    dis = distance3d(n1, n2)
    return StaticEdge(src, des, id, g, dis)
end

function update_src_des(e::StaticEdge, src::Int, des::Int)
    return StaticEdge(src, des, e.id, e.grade, e.distance)
end

function Base.getproperty(e::StaticEdge, name::Symbol)
    if name == :grades
        return (e.grade,)
    elseif name == :distances
        return (e.distance,)
    elseif name == :nodes
        return SA[e.src, e.des]
    else
        return getfield(e, name)
    end
end

"""
A way is a road with multiple road segments (encodes in nodes)
"""
@with_kw mutable struct Way <: AbstractWay
	src::Int64
	des::Int64
	id::Int64
	islink::Bool = false # if this way is a link that connect two roads
	nodes::Vector{Int64}
    grades::Vector{Float64}
    distances::Vector{Float64}
end

islink(way::Way) = way.islink

function merge(e1::Way, e2::Way)
	error("unimplemented")	
end

function Way(net, src::Int64, des::Int64, id::Int64, islink::Bool, nodes::Vector{Int64})
    node_v = getnode.((net,), nodes)
    grade_v = grades(node_v)
    dis_v = distance3d.(node_v[1:end-1], node_v[2:end])
    return Way(src, des, id, islink, nodes, grade_v, dis_v)
end