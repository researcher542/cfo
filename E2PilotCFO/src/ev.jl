
"""
The electric vehicle's model to caculate the energy consumption.
It is ok to keep it mutable since we are no much of them...
"""
@with_kw mutable struct EV{VT, opi_T,msi_T,sg_T} <: ep.AbstractVehicle
	m::Float64 = 36000. # mass in kg
	Af::Float64 = 10. # front area in m^2
	Cd::Float64 = 0.36 # drag coefficient
	pwr_out_perc::VT = SA[0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.20, 0.40, 0.60, 0.80, 1.00] # 
	eff_map::VT = SA[0.83, 0.85, 0.87, 0.89, 0.90, 0.91, 0.93, 0.94, 0.94, 0.93, 0.92]  # efficient map
	Cr::Float64 = 0.0094 # the rolling coefficient
	eta::Float64 = 0.85 # driveline efficiency, 0.94 in the paper, but 0.85 used in FastSim
	aux_pow::Float64 = 10.0 * 1e3 # auxliary power, in W, accounts for the power when idling
	max_power::Float64 = 380.0 * 1e3 # max power, in W
    cap::Float64 = ep.kwh2j(1000.0) #  battery capacity in joule
    regen_eff::Float64 = 1/0.98 # The efficiency of regenerative system. It should be greater than 1 to make the conculation consistent.
    objtype::AbstractObjective = ObjEnergy()

    # computed properties
    out_power_itp::opi_T = missing # store the information of output power in sitp
    max_spd_itp::msi_T = missing # store the information of maximum speed subject to maximum power limit
    max_spd_vec::Union{Vector{Float64}, Missing} = missing # store the information in a data
	sol_grid::sg_T = missing
end

function Base.show(io::IO, ::Type{T}) where T <: EV
    Base.show(io, "EV")
end

function is_no_regen(ev::EV)
    return isinf(ev.regen_eff) 
end

function set_no_regen!(ev::EV)
    ev.regen_eff = Inf 
    ev = ep.update_veh!(ev)
    return ev
end

function get_ev()
    kwargs = ep.compute_prop_veh(EV)
    return EV(;kwargs...)
end

function ep.energy_cost_acc(veh::EV, nodes::Vector{S}, spd_v::Vector{Float64})::Float64 where S<: ep.AbstractNode
    return 0.0
end

# """
# dynamically compute and load the solution grid
# """
# function Base.getproperty(veh::EV, sym::Symbol)
#     #if sym == :out_power_itp && ismissing(getfield(veh, sym))
#     #    veh.out_power_itp = get_out_power_itp(veh)
#     if (sym == :max_spd_itp || sym == :max_spd_vec) && ismissing(getfield(veh, sym))
#         update_spd_limit_data!(veh)
#     end
#     if sym == :sol_grid && ismissing(getfield(veh, sym))
#         update_solgrid!(veh)
#     end
#     return getfield(veh, sym)
# end

"""
required power to maintain the constant speed
return: power in W
"""
function ep.required_power(ev::EV, theta::Real, spd::Real, H::Real=0.0)
	power = ep.resistance_force(ev, theta, spd, H)*spd
	return power 
end

ep.efficiency(ev::EV, power::Float64) = efficiency_ev(ev, power)

"""
return the efficiency of the vehicle given at certain power
"""
function efficiency_ev(ev::EV, power::Float64)
    if power < 0.0
        return ev.regen_eff
    else
        return ep.efficiency_raw(ev, power)
    end
end

ep.augmentedtime(net::Network, src::Int, des::Int, lam::Float64, ev::EV) = ep.augmentedtime_raw(net, src, des, lam, ev)


"""
The power of fuel consumption

return: in W
"""
function ep.output_power_raw(ev::EV, theta::Float64, spd::Float64, H=0.0; nan_flag=false)
	pow_out = ep.required_power(ev, theta, spd, H)
	power = pow_out / ev.eta + ev.aux_pow
	if power > ev.max_power
		if nan_flag
			return Inf
		end
	end
	eff = efficiency_ev(ev, power)
	pow = power / eff
	return pow
end

ep.output_power(veh::EV, θ::Float64, spd::Float64, H::Float64=0.0) = ep.output_power_raw(veh, θ, spd, H)

function get_energy_cost(veh::EV,  way::Way, spd::Float64, H::Float64 = 0.0)
    dis::Float64 = way.distances[1]
    theta::Float64 = way.grades[1]
    pow::Float64 = output_power(veh, theta, spd, H)
    t::Float64 = dis / spd 
    return pow * t
end

get_energy_cost_inv(net::Network, ev::EV, src::Int, des::Int, cost0::Float64) =  get_energy_cost_inv(net::Network, ev::EV, getway(net, src, des), cost0::Float64)

"""
Given cost, return the travelling time
"""
function get_energy_cost_inv(net::Network, ev::EV, way::Way, cost0::Float64)
    dis::Float64 = way.distances[1]
    (;src, des) = way
    min_spd, max_spd = get_minmax_speed(net, src, des, ev)
    function f(spd)
        energy_cost = get_energy_cost(ev, way, spd)
        return energy_cost - cost0
    end

    method = SA[Roots.A42(), Roots.Bisection()][1]
    opt_spd = find_zero(f, (min_spd, max_spd), method; xrtol=1e-4)
    e_cost = get_energy_cost(ev, way, opt_spd)
    @assert ( abs(e_cost) <= abs(cost0) * (1+ 1e-4) + 1e-6)  
    t = dis / opt_spd
    return t
end