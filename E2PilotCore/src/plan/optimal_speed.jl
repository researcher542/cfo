
optimal_speed(veh::AbstractVehicle, theta::Real,
	lambda::Real, min_spd::Real=g_min_highway_speed,
	max_spd::Real=g_max_highway_speed) = optimal_speed(veh.sol_grid, theta, lambda, min_spd, max_spd)

function optimal_speed(sol_grid::SolutionGrid{T}, theta::Float64,
	lambda::Real, min_spd::Real=g_min_highway_speed,
	max_spd::Real=g_max_highway_speed) where T

	# theta_idx = searchsortedlast(sol_grid.theta_range, theta)
	theta_idx = search_over_range(sol_grid.theta_range, theta)
	if theta_idx == length(sol_grid.theta_range) + 1
		#@warn "theta=$theta exceed the maximum grade"
		theta_idx -= 1
	elseif theta_idx == 0
		theta_idx += 1
		#@warn "theta=$theta exceed the minimum grade"
	end
    
    # it seems tricky, extrapolation seems better than Interpolations
    # itp = sol_grid.itp_v[theta_idx]
    
    lam0::Float64 = lambda
    # itp = sol_grid.itp_v[theta_idx]
    # opt_spd::Float64 = itp(lam0)
    lam_vec = @view sol_grid.sol_mat[:,theta_idx]
    opt_spd::Float64 = -1.0
    spd_range = sol_grid.spd_range
    idx = searchsortedlast(lam_vec, lam0)
    opt_spd = spd_range[idx]
    # @inbounds for ispd in length(spd_range):-1:1
    #     lam = lam_vec[ispd] 
    #     if lam <= lam0
    #         opt_spd = spd_range[ispd]
    #     end
    # end

    
    opt_spd = clamp(opt_spd, min_spd, g_max_highway_speed)
    # maximum speed is stricter than min_spd
    opt_spd = clamp(opt_spd, 0.0, max_spd)
    
    if opt_spd == 0
        @warn "spd==0"
        # @infiltrate
    end

	return opt_spd
end

"""
search over a range with uniform steps, and increasing order
"""
function search_over_range(range, x::Float64)::Int
    n = length(range)
    min = range[1]
    max = range[end]
    x = clamp(x, min, max) 
    ## normalize x to [0,1]
    normal_x = (x-min) / (max-min)
    return round(Int, normal_x*n)
end

