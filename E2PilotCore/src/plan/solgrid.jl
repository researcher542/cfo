"""
"""


@with_kw struct SolutionGrid{T}
	spd_range::Vector{Float64}
	theta_range::Vector{Float64}
	sol_mat::Matrix{Float64}
	itp_v::T # Interpolations object
	max_lam_vec::Vector{Float64}
end

getrange(itp::Interpolations.Extrapolation) = itp.itp.knots[1]

# function Base.show(io::IO, sg::SolutionGrid)
#     Base.print(io, "Solution grid with spd_range $(sg.spd_range).\n theta_range $(sg.theta_range)")
# end

SolutionGrid(veh::AbstractVehicle; kwargs...) = SolutionGrid(veh, collect(g_theta_range), collect(g_spd_range) ; kwargs...)

function max_lambda(sol_mat)
	f(x::Float64)::Float64 = isinf(x) ? -Inf64 : x
	max_lam = maximum(f, sol_mat)
	return max_lam
end
max_lambda(sol_grid::SolutionGrid) = maximum(sol_grid.max_lam_vec)
max_lambda(veh::AbstractVehicle) = max_lambda(veh.sol_grid)

function SolutionGrid(veh::AbstractVehicle, theta_range::Vector{Float64}, spd_range::Vector{Float64}; kwargs...)
	sol_mat = solution_matrix(veh, theta_range, spd_range; objtype=veh.objtype, kwargs...)
    itp_v = get_itp_v(sol_mat, theta_range, spd_range)
	max_lam_vec = [max_lambda(sol_mat[:,idx]) for idx in 1:size(sol_mat,2)]
	return SolutionGrid(spd_range, theta_range, sol_mat, itp_v, max_lam_vec)	
end

"""
g(t) = t*f(D/t) + λ t, where f(v) is the fuel rate function in W. We want the optimal speed v such that:
g'(t) = f(v) - v*f'(v) + λ == 0.
Let h(v) = v*f'(v) - f(v) == λ, 
we can build a grid (v,θ) that contain the value of h(θ,v). Latter, given λ, we can retrive the value of v*
"""
function solution_matrix(veh::AbstractVehicle, theta_range = g_theta_range, spd_range = g_spd_range; objtype::AbstractObjective = ObjEnergy(), kwargs...)
	mat = Matrix{Float64}(undef ,length(spd_range), length(theta_range))
    for (itheta,theta) in enumerate(theta_range)
		for (ispd,spd) in enumerate(spd_range)
            h = spd2lambda(veh, theta, spd, objtype; kwargs...)
			if h > 0
				#@debug spd, theta, fd, h
			end
			mat[ispd, itheta] = isnan(h) ? Inf : h
		end
	end
	return mat
end

"""
The map from speed to lambda. See the document of solution_matrix
"""
function spd2lambda(veh::AbstractVehicle, theta::Real, spd::Real, objtype::AbstractObjective; kwargs...)
    fval = objective(veh, theta, spd, objtype; kwargs...)::Float64
    # derivative of the objective
	fd = objective_d(veh, theta, spd, objtype; kwargs...)::Float64
	h = spd*fd - fval
    return h
end

"""
For Emission objective, the map from speed to lambda is a little tricky.
Because the emission function is not continuous.
Note: but the emission objective is continuous.
"""
function spd2lambda(veh::AbstractVehicle, theta::Real, spd::Real, objtype::ObjEmission; kwargs...)
    kw_dict = Dict(kwargs)
    swi_spd_vec = haskey(kw_dict, :swi_spd_vec) ? kw_dict[:swi_spd_vec] : g_switching_speed_vec
    args = [veh, theta, spd, objtype]
    fval = objective(args...; kwargs...)
    fd = objective_d(args...; kwargs...)
    #rpm = speed2rpm(spd)
    #eps = 1e-8
    #for swi_spd in swi_spd_vec
    #    # if the speed is near the swiching speed, should be careful with the derivative
    #    if abs(rpm-swi_spd) <= eps
    #        spd_eps = (rpm <= swi_spd) ? spd - eps : spd + eps
    #        args[3] = spd_eps
    #        fval_eps = objective(args...; kwargs...)
    #        fd = (fval-fval_eps)/eps
    #        fd = fd * ((rpm <= swi_spd) ? 1 : -1)
    #    end
    #end
	h = spd*fd - fval
    return h
end


"""
Given the solution matrix that stores mat[spd,theta] = h(v,theta)
Here, we deine an interpolator such that
h_theta^-1(λ) = v^star
Note that due to the jump property of the fuel rate function, same λ may have multiple v^star, here we pick the maximum one
"""
function get_itp_v(sol_mat::Matrix, theta_range, spd_range)
    itp_v = Vector{Interpolations.Extrapolation}(undef, length(theta_range))
    
    for (itheta, theta) in enumerate(theta_range)
        lam_vec = @view sol_mat[:,itheta]
        lam_vec_new = Vector{Float64}()
        spd_vec_new = Vector{Float64}()
        n_lam = length(lam_vec)
        # 
        for lam in lam_vec
            if isinf(lam) | (lam < 0) continue end
            # search from the end
            for idx in n_lam:-1:1
                flag1 = lam >= lam_vec[idx]
                # flag2 = !(lam in lam_vec_new)
                if flag1 
                    # @debug lam
                    spd_ = spd_range[idx]
                    push!(lam_vec_new, lam)
                    push!(spd_vec_new, spd_)
                    break
                end
            end
        end
        if length(lam_vec_new) == 1
            lam_vec_new = [lam_vec_new[1], lam_vec_new[1]+1]
            spd_vec_new = [spd_vec_new[1], [spd_vec_new[1]]]
            @debug lam_vec_new, spd_vec_new
        elseif isempty(lam_vec_new) 
            max_spd = maximum(spd_range)
            lam_vec_new = [0, 1]
            spd_vec_new = [max_spd, max_spd]
        end
        I = sortperm(lam_vec_new)
        # @debug "$theta, $(isempty(lam_vec_new))"
        # itp_v[itheta] = LinearInterpolation(lam_vec_new[I], spd_vec_new[I], extrapolation_bc=Line())
        lam_vec_new1 = lam_vec_new[I]
        spd_vec_new1 = spd_vec_new[I]
        Interpolations.deduplicate_knots!(lam_vec_new1)
        itp_v[itheta] = linear_interpolation(lam_vec_new1, spd_vec_new1, extrapolation_bc=Line())
    end # end of theta loop
    # make it the concrete type 
    return [itp for itp in itp_v]
end

# Wrappers for different objective types

"""
The fuel rate function in kW
"""
objective(veh::AbstractVehicle, theta::Real, speed::Real, objtype::ObjEnergy; kwargs...) = output_power(veh, theta, speed)
# objective_d(veh::AbstractVehicle, theta::Real, speed::Real, objtype::ObjEnergy) = fuel_power_d(veh, theta, speed)
function objective_d(veh::AbstractVehicle, theta::Real, speed::Real, objtype::AbstractObjective; kwargs...) 
    eps	= 1e-10
	obj1 = objective(veh, theta, speed-eps, objtype; kwargs...)	
	obj2 = objective(veh, theta, speed+eps, objtype; kwargs...)	
	return (obj2 - obj1) / (2*eps)
end