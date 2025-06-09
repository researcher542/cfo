
abstract type AbstractNode end


"""
A immutable struct for the node data. Basically
"""
@with_kw struct StaticNode <: AbstractNode
	lat::Float64
	lon::Float64
	ele::Float64	# elevation
	id::Int64
	idx::Int64 	# The idx in the graph
    name::Symbol = Symbol()
end

function update_idx(nd::StaticNode, idx::Int64)
    new_nd = StaticNode(nd.lat, nd.lon, nd.ele, nd.id, idx, nd.name)
    return new_nd
end

@with_kw mutable struct Node <: AbstractNode
	lat::Float64
	lon::Float64
	ele::Float64	# elevation
	id::Int64
	idx::Int64 	# The idx in the graph
    name::Symbol = Symbol()
	function Node(lat::Union{String, Float64}, lon::Union{String, Float64}, ele::Union{String, Float64}, id::Int64, idx::Int64, name::Symbol = Symbol())
		f = x -> (isa(x, String) ? parse(Float64, x) : x)
		return new(f(lat), f(lon), f(ele), id, idx, name)
	end
end
# Base.getproperty(node::AbstractNode, sym::Symbol) = getfield(node, sym)

isless(n1::LatLon, n2::LatLon) = (n1.lat < n2.lat) || (n1.lat == n2.lat && (n1.lon < n2.lon))

Node(lat::Real, lon::Real, ele::Real, id::Int64, idx::Int64) = Node(promote(lat,lon,ele)..., id, idx)
LLA(n::AbstractNode) = LLA{Float64}(n.lat, n.lon, n.ele) 
LatLon(n::AbstractNode) = LatLon{Float64}(n.lat, n.lon)

const haversine = Haversine()
distance2d(n1, n2) = haversine((n1.lon, n1.lat), (n2.lon, n2.lat)) 
distance3d(n1::T, n2::T) where T <: Union{AbstractNode, LLA} = sqrt( (n1.ele - n2.ele)^2 + distance2d(n1, n2)^2 )

#distance2d(n1::LatLon, n2::LatLon) = euclidean_distance(n1, n2)
#distance2d(n1::AbstractNode, n2::AbstractNode) = euclidean_distance(LatLon(n1), LatLon(n2))
#distance3d(n1::AbstractNode, n2::AbstractNode) = euclidean_distance(LLA(n1), LLA(n2))
distance(nodes::Vector{S}) where S <: AbstractNode = sum(distance3d.(nodes[1:end-1], nodes[2:end]))

"""
The grade/slope between to coordinates
return: unitless ratio of tan(theta)
"""
function grade(n1::AbstractNode, n2::AbstractNode)::Float64
	dis2d = distance2d(n1,n2)
	# if it is less than 1 meter
	if dis2d <= 1
		return 0.0
	end
	del_ele = n2.ele - n1.ele
    g = del_ele/dis2d
	if abs(g) > 100
		@warn "irrugular grade $g"
		@show del_ele
		@show dis2d
	end
	return g
end

function grades(nodes::Vector{S})::Vector{Float64} where S<:AbstractNode
	@assert length(nodes) >= 2
    n = length(nodes)
    grade_v = zeros(Float64, n-1)
    for i in 1:n-1
        grade_v[i] = grade(nodes[i], nodes[i+1])
    end
    return grade_v
end

function polyline_src(polyline::String)::LatLon
	pts = decodePolyline(polyline)
	return LatLon(pts[1,1], pts[1,2])
end

function polyline_des(polyline::String)::LatLon
	pts = decodePolyline(polyline)
	return LatLon(pts[end,1], pts[end,2])
end

function pts2latlon(pts::Matrix{T})::Vector{LatLon{T}} where T
	latlons = LatLon.(pts[:,1], pts[:,2])
	return latlons
end

function latlon2pts(latlons::Vector)::Matrix{Float64} 
	#@debug "latlon2pts called"
	n = length(latlons)
	pts = zeros(n, 2)
	for i in 1:n
		latlon = latlons[i]
		pts[i,1] = latlon.lat
		pts[i,2] = latlon.lon
	end
	return pts
end
nodes2pts = latlon2pts

# function reduce_by_distance(latlons::Union{Vector{LatLon{T}},Vector{Node}}; min_dis=k_min_nodes_distance) where T <: Real
function reduce_by_distance(latlons; min_dis=k_min_nodes_distance) 
	if length(latlons) <= 2
		return latlons
	end
	latlons_red = [latlons[1]]
	for latlon in latlons[2:end-1]
		if distance2d(latlon, latlons_red[end]) > min_dis
			push!(latlons_red, latlon)
		end
	end
	# push!(latlons_red, latlons[end])
	latlons_red = [latlons_red[1], latlons_red[2:end-1]..., latlons[end]]
	return latlons_red
end

"""
"""
function reduce_by_distance(pts::Matrix; kwargs...)::Matrix
	if size(pts,1) <= 2
		return pts
	end
	latlons = pts2latlon(pts)
	latlons_red = reduce_by_distance(latlons;kwargs...)
	pts_red = latlon2pts(latlons_red)
	return pts_red
end

function is_valid_grade(grade::Float64)
	return grade >= minimum(g_theta_range) && grade <= maximum(g_theta_range)
end

"""
nodes: vector of nodes with grade implemented
"""
function reduce_by_grade(nodes::AbstractArray)
	n_node = length(nodes)
	if length(nodes) <= 2
		return nodes
	end
	tol = (g_theta_range[2] - g_theta_range[1])
	nodes_red = [nodes[1]]
	cur_grade = 1e3
	for (node1,node2) in zip(nodes[1:end-1],nodes[2:end-1])
		tmp_grade = grade(node1, node2)
		if !is_valid_grade(tmp_grade) 
			continue
		end
		if abs(cur_grade - tmp_grade) > tol 
			push!(nodes_red, node2)
			cur_grade = tmp_grade
		end
	end
	push!(nodes_red, nodes[end])
	# @debug "reduce nodes from $n_node to $(length(nodes_red)) by grade"
	return nodes_red
end


"""
Given a vector of nodes, check if their grades are feasible
"""
function check_grades(nodes::Vector{Node})
	max_theta = maximum(g_theta_range)
	min_theta = minimum(g_theta_range)
	grade_v = grade.(nodes[1:end-1], nodes[2:end])
	f = theta -> (abs(theta) > max_theta)
	index = findall(f, grade_v)
	dis_v = distance2d.(nodes[index], nodes[index.+1])
	grade1_v = grade_v[index]
	@debug index	
	@debug dis_v
	@debug grade1_v
	#@debug grade_v
end

function plot_ele(nodes::Vector{Node}, data_v=missing)
	dis_v1 = distance2d.(nodes[1:end-1], nodes[2:end])
	dis_v = accumulate(+, [0; dis_v1])
	if ismissing(data_v)
		ele_v = [n.ele for n in nodes]
	else
		ele_v = data_v
	end
	plot(dis_v, ele_v)
end



