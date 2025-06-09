#=
The module that related to vehicle and fuel model information
=#

const k_rho_a = 1.2256  # air density, in kg/m^3
const k_gravity_acc = 9.8066 # gravity acceleration

"""
The vehicle's model to caculate the fuel consumption.
refer to: Fuel consumption model for heavy duty diesel trucks: Model development and testing.
It is ok to keep it mutable since we are no much of them...
"""
@with_kw mutable struct Vehicle{TA, TB, TC, TD} <: AbstractVehicle
	m::Float64 = 36000. # mass in kg
	Af::Float64 = 8.5502 # front area in m^2
	Cd::Float64 = 0.7 # drag coefficient
	pwr_out_perc::TA = SA[0, 0.005, 0.015, 0.04, 0.06, 0.10, 0.14, 0.20, 0.40, 0.60, 0.80, 1.00] # 
	eff_map::TA = SA[0.10, 0.14, 0.20, 0.26, 0.32, 0.39, 0.41, 0.42, 0.41, 0.38, 0.36, 0.34] # efficient map
	Cr::Float64 = 0.0094 # the rolling coefficient
	eta::Float64 = 0.85 # driveline efficiency, 0.94 in the paper, but 0.85 used in FastSim
	aux_pow::Float64 = 10.0 * 1e3 # auxliary power, in w, accounts for  the power when idling
	max_power::Float64 = 380.0 * 1e3 # max power, in W
    objtype::AbstractObjective = ObjEnergy()
    
    # computed properties
    out_power_itp::TB = missing # store the information of output power in sitp
    max_spd_itp::TC = missing # store the information of maximum speed subject to maximum power limit
    max_spd_vec::Union{Vector{Float64}, Missing} = missing # store the information in a data
	sol_grid::TD = missing

end


"""
cons: the constructor
"""
function compute_prop_veh(cons)
    veh = cons()
    kwargs = Dict()
    # @debug "getting out_power_itp"
    kwargs[:out_power_itp] = get_out_power_itp(veh)
    veh = cons(;kwargs...)
    # @debug "getting max_spd itp"
    max_spd_itp, max_spd_vec = get_spd_limit_data(veh)
    kwargs[:max_spd_itp] = max_spd_itp
    kwargs[:max_spd_vec] = max_spd_vec
    veh = cons(;kwargs...)
    sol_grid = SolutionGrid(veh)
    # @debug "getting sol grid"
    kwargs[:sol_grid] = sol_grid
    # @debug "sol grid got."
    # :sol_grid => sol_grid
    # veh = typeof()
    return kwargs
end

function get_veh()
    kwargs = compute_prop_veh(Vehicle)
    return Vehicle(;kwargs...)
end

function get_spd_limit_data(veh::AbstractVehicle)
    #theta_range = g_min_theta:1e-4:g_max_theta
    theta_range = g_theta_range
    max_spd_vec = [ get_max_speed_raw(veh, θ)  for θ in theta_range]
    #itp = LinearInterpolation(theta_range, max_spd_vec, extrapolation_bc=Line()) 
    itp = LinearInterpolation(theta_range, max_spd_vec) 
    #veh.max_spd_itp = itp
    #veh.max_spd_vec = max_spd_vec
    return itp, max_spd_vec
end

function get_max_speed(veh::AbstractVehicle, θ::Float64)
    θ::Float64 = clamp(θ, g_min_theta, g_max_theta)
    itp = veh.max_spd_itp
    return itp(θ)
end

is_undef_pow(pow::Float64) = (isinf(pow) || isnan(pow))

"""
get maximum speed according to maximum power
"""
function get_max_speed_raw(veh::AbstractVehicle, θ::Real, H=0.0)
    fmax = output_power(veh, θ, g_max_highway_speed, H)
    if !is_undef_pow(fmax) 
        return g_max_highway_speed
    end
    min_spd = 1e-3
    fmin = output_power(veh, θ, min_spd, H)
    if is_undef_pow(fmin)
        return min_spd
    end
    function func(spd::Float64)
        pow = output_power(veh, θ, spd)
        return is_undef_pow(pow)
    end
    #spd_range = g_min_spd:0.1:g_max_spd
    spd_range = g_spd_range
    idx = findfirst(func, spd_range) 
    # idx should minus 2 to avoid interpolate between a valid number and inf power
    idx1 = (isnothing(idx) || idx <= 2) ? 1 : idx - 2
    return spd_range[idx1]
end

"""
A combination of traffic limit and grade speed limit
"""
function get_max_speed(net::Network, n1::Int64, n2::Int64, veh::AbstractVehicle)
    θ = grade(net, n1, n2)
    spd1::Float64 = get_max_speed(net, n1, n2)
    spd2::Float64 = get_max_speed(veh, θ)
    spd = ifelse(spd1 < spd2, spd1, spd2)
    return spd
end

function get_min_speed(net::Network, n1::Int, n2::Int, veh::AbstractVehicle)
    return g_min_highway_speed 
end

function get_minmax_speed(net::Network, n1::Int, n2::Int, veh::AbstractVehicle)
    max_spd = get_max_speed(net, n1, n2, veh)
    min_spd = get_min_speed(net, n1, n2, veh)
    min_spd = clamp(min_spd, g_min_spd, max_spd)
    return min_spd, max_spd
end

function get_minmax_t(net::Network, u::Int, v::Int, veh::AbstractVehicle)
    min_spd, max_spd = get_minmax_speed(net, u, v, veh) 
    dis = distance(net, u, v)
    return (dis/max_spd, dis/min_spd)
end

# update_solgrid!(veh::AbstractVehicle; kwargs...) = 0.0
function update_solgrid!(veh::AbstractVehicle; kwargs...)
    sol_grid = SolutionGrid(veh; kwargs...)
    veh.sol_grid = sol_grid
    return veh
end

function update_veh!(veh::AbstractVehicle; kwargs...)
    update_solgrid!(veh; kwargs...)
    veh.out_power_itp = get_out_power_itp(veh)
    max_spd_itp, max_spd_vec = get_spd_limit_data(veh)
    veh.max_spd_itp = max_spd_itp
    veh.max_spd_vec = max_spd_vec
    return veh
end

function set_weight!(veh::Vehicle, weight::Float64)
    if veh.m == weight
        return veh
    end
    veh.weight = weight
    update_veh!(veh)
    return veh
end

function set_objtype!(veh::AbstractVehicle, objtype::AbstractObjective; kwargs...)
    # @debug "setting veh obj $objtype..."
    veh.objtype = objtype
    update_veh!(veh; kwargs...)
    return veh
end

set_objtype(veh::Vehicle, objtype::AbstractObjective; kwargs...) = set_objtype!(deepcopy(veh), objtype; kwargs...)

"""
Get the maximum energy cost perge edge for the network.
"""
function get_max_e_cost(net::AbsNet, veh::AbstractVehicle)
    max_cost = 0.0
    iway = -1
    for (k, w) in net.waydata
        u = w.src
        v = w.des
        (min_t, max_t) = get_minmax_t(net, u, v, veh)
        t = min_t
        e_cost = energy_cost_on_road(net, u, v, veh, t)
        if max_cost < e_cost
            max_cost = e_cost
            iway = k
        end
    end
    return max_cost
end


include("emission.jl")
include("energy.jl")

#const k_veh_dict = get_veh_dict()
#const k_default_veh = get_vehicle()