"""
calculate the energy summary
"""

function energy_cost_on_road(net::Network, u::Int, v::Int, veh::AbstractVehicle, t::Float64)
    (u == v) && return 0.0
    dis = distance(net, u, v)
    θ = grade(net, u, v)
    spd = dis / t
    pow = output_power(veh, θ, spd)
    cost = pow * t
	return cost::Float64
end

"""
caculate the energy consumption for a path and speed profile
	return: energy in Joule
"""
function energy_cost_on_road(veh::AbstractVehicle, nodes::Vector{S}, spd_v::Vector{Float64}, nan_flag::Bool=false)::Float64  where S <: AbstractNode
	dis_v = distance3d.(nodes[1:end-1], nodes[2:end])
	grade_v = grades(nodes)
	div(a,b) = (b == 0) ? 0.0 : a/b
	time_v = div.(dis_v, spd_v)
	# f(g, spd) = output_power(veh, g, spd)
    n = length(dis_v)
    #power_v = Vector{Float64}(undef, )
    cost::Float64 = 0.0
    for i in 1:n
        pow = output_power_raw(veh, grade_v[i], spd_v[i]; nan_flag = nan_flag)
        # if !nan_flag
        # else
        #     pow = output_power(veh, grade_v[i], spd_v[i])
        # end
        # @assert !isinf(pow)
        cost += time_v[i] * pow
    end
	return cost
end

"""
Energy cost by acceleration/changing speed
"""
function energy_cost_acc(veh::AbstractVehicle, theta::Float64, spd1::Float64, spd2::Float64)::Float64
	avg_spd = (spd1 + spd2) / 2
	power = output_power_raw(veh, theta, avg_spd)
	eta = efficiency(veh, power)
	e =  1/2.0 * veh.m * abs(spd1^2 - spd2^2) / eta
	return e
end

function energy_cost_acc(veh::AbstractVehicle, nodes::Vector{S}, spd_v::Vector{Float64})::Float64 where S<:AbstractNode
    if length(nodes) <= 2
        return 0.0
    end
	grade_v = grades(nodes)
	energy_cost_v = energy_cost_acc.((veh,), grade_v[1:end-1], spd_v[1:end-1], spd_v[2:end] )
	return sum(energy_cost_v)/2
end

function total_energy_cost(veh::AbstractVehicle, nodes::Vector{S}, spd_v::Vector{Float64}, nan_flag::Bool = false) where S <: AbstractNode
	e1 = energy_cost_on_road(veh, nodes, spd_v, nan_flag)	
	e2 = energy_cost_acc(veh, nodes, spd_v)
    e2 = is_undef_pow(e2) ? 0.0 : e2
	if e2 != 0
		#@debug "road to speedchange ratio = $(e1/e2)" 
	end
	cost = e1 + e2
	return cost
end

function cost(objtype::AbstractObjective, veh::AbstractVehicle, net::Network, i::Int, j::Int, spd::Float64; kwargs...)
    nodes = [getnode(net, i), getnode(net, j)]
    return total_cost(veh, nodes, [spd]; kwargs...) 
end

function cost(objtype::ObjEnergy, veh::AbstractVehicle, net::Network, i::Int, j::Int, spd::Float64)
    theta = grade(net, i, j) 
    dis::Float64 = distance(net, i, j)
    t::Float64 = dis/spd
    pow = output_power(veh, theta, spd)
    return pow*t
end

"""
high-level summary of cost, can
"""
function total_cost(veh::AbstractVehicle, nodes::Vector{S}, spd_v::Vector{Float64}; param_vec=k_emission_param_vec, swi_spd_vec=g_switching_speed_vec) where S <: AbstractNode
    objtype = veh.objtype
    if isa(objtype, ObjEnergy)
        return total_energy_cost(veh, nodes, spd_v)
    elseif isa(objtype, ObjEmission)
        return total_emission_cost(veh, nodes, spd_v, param_vec, swi_spd_vec)
    else 
        @warn "cost summary for $(typeof(objtype)) not implemented"
        return 0.0
    end
end


