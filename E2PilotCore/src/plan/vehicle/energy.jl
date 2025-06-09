
using Polynomials


"""
fuel cost related functions
"""

"""
theta: unitless theta = diff_ele/ plane_distance

spd: in m/s

H: the altitude in km

Refer to the source code of the FastSim 

return: Force in N
"""
function resistance_force(v::AbstractVehicle, theta::Float64, spd::Float64, H::Real=0.0)
	# the correction factor for height, 
	# should be 1 - 0.085 H, H is the altitude in km
	Ch = 1 - 0.085 * H
	f_air = 1/2 * k_rho_a * v.Cd * Ch * v.Af * spd^2
	f_roll = v.m * k_gravity_acc * v.Cr
	f_gravity = v.m * k_gravity_acc * sin(atan(theta))
	F = f_air + f_roll + f_gravity
	return F
end

"""
required power to maintain the constant speed
return: power in W
"""
function required_power(v::Vehicle, theta::Real, spd::Real, H::Real=0.0)
	power = resistance_force(v, theta, spd, H)*spd
    power = power < 0 ? 0.0 : power
	return power 
end

get_eff_map(veh::AbstractVehicle) = veh.eff_map
get_pwr_out_perc(veh::AbstractVehicle) = veh.pwr_out_perc



efficiency(veh::AbstractVehicle, power::Float64) = efficiency_raw(veh, power)

"""
return the efficiency of the vehicle given at certain power
"""
function efficiency_raw(veh::AbstractVehicle, power::Float64)
    (;max_power::Float64) = veh
    # Function barrier to enforce type stability.
    eff_map = get_eff_map(veh)
    pwr_out_perc = get_pwr_out_perc(veh)
	perc::Float64 = power / max_power
	if perc > 1.0
		#@warn "required power $(power) larger that maxpower $(v.max_power)"
		perc = 1.0
	end
	#vec = [x >= percent for x in v.pwr_out_perc]
	#idx = argmax(vec) - 1
	#idx = (idx == 0) ? 1 : idx

    # idx_tmp1 = findfirst(x-> x>=perc, pwr_out_perc) 
    # idx_tmp::Int = (isnothing(idx_tmp1) || idx_tmp1 == 1) ? 1 : idx_tmp1 - 1
    # idx::Int = idx_tmp
    @assert (power >= 0)
    ## self implemented findfirst for performance
    idx::Int = -1
    for iperc in 1:length(pwr_out_perc)-1
        perc1 = pwr_out_perc[iperc]
        perc2 = pwr_out_perc[iperc+1]
        if perc1 <= perc <= perc2
            idx = iperc
            break
        end
    end
    
  
	y1::Float64 = eff_map[idx]
	y2::Float64 = eff_map[idx+1]
	x1::Float64 = pwr_out_perc[idx]
	x2::Float64 = pwr_out_perc[idx+1]
	eff = y1 + (y2-y1)/(x2-x1)*(perc-x1)
	return eff
end

"""
The power of fuel consumption

return: in W
"""
function output_power_raw(v::AbstractVehicle, theta::Float64, spd::Float64, H=0.0; nan_flag=true)
	pow_out = required_power(v, theta, spd, H)
	power = pow_out / v.eta + v.aux_pow
	if power > v.max_power
		if nan_flag
			return Inf
		end
	end
	eff = efficiency(v, power)
	pow = power / eff
	return pow
end

"""
Note: this inexact solution might cause negative cycle for EV
"""
function output_power(veh::AbstractVehicle, θ::Float64, spd::Float64, H=0.0)  
    θ0::Float64 = clamp(θ, g_min_theta, g_max_theta)
    spd0::Float64 = clamp(spd, g_min_spd, g_max_spd)
    @inbounds pow::Float64 = veh.out_power_itp(θ0, spd0)
    # pow = output_power_raw(veh, θ0, spd0, H)
    return pow
end

function get_out_power_itp(veh::AbstractVehicle)
    spd_range = g_min_spd:0.01:g_max_spd
    theta_range = g_min_theta:1e-4:g_max_theta
    A = [output_power_raw(veh, θ, spd) for θ in theta_range, spd in spd_range]
    itp = interpolate(A, (BSpline(Linear()), BSpline(Linear())))
    sitp = scale(itp, theta_range, spd_range)
    return sitp
end

"""
Derivative of fuel power function f(v) with respective to speed v.
	Numerically computed.
"""
function output_power_d(v::AbstractVehicle, theta::Float64, spd::Float64, H=0.0)
	eps	= 1e-8
	pow1 = output_power(v, theta, spd-eps, H)	
	pow2 = output_power(v, theta, spd+eps, H)	
	return (pow2 - pow1) / (2*eps)
end

"""
For now, just simply simluate for all the speeds and perform a regression

theta: unitless theta = diff_ele/plane_distance

H: the altitude in km

return: A Polynomials obj

	f(v) = Poly(v)
	unit: kW
"""
function fuel_func(v::Vehicle, theta::Float64, H::Float64=0.0; plot_flag=false)
	spds = copy(g_spd_range)
	f = spd -> output_power(v, theta, spd, H)

	fuel_powers = f.(spds)
	## some of the spd may exceed the maximum power, we may disgrad those data
	idx = findfirst(isnan, fuel_powers)
	max_idx = isnothing(idx) ? length(spds) : idx -1
	spds = spds[1:max_idx]
	fuel_powers = fuel_powers[1:max_idx]

	# For better accuracy.
	tol = 1e-2
	deg = 3
	fuel_poly = Polynomial()
	fit_y = 0
	while deg < length(spds) - 1
		fuel_poly = fit(spds, fuel_powers, deg)
		fit_y = fuel_poly.(spds)
		err = fit_y - fuel_powers
		rel_err = err./ fuel_powers
		rel_err = LA.norm(rel_err, 2) / sqrt(length(rel_err))
		if rel_err < tol
			break
		else
			#@debug "for theta=$(theta*100)% deg=$deg relative error is $(rel_err*100)%"
		end
		deg += 1
	end

	if plot_flag
		p = plot()
		plot!(p, spds, fuel_powers, label="original data")
		#plot!(p, fuel_poly,label="poly")
		plot!(p, spds, fit_y, label="poly")
		display(p)
	end
	return fuel_poly
end

