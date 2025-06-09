"""
"""

"""
Refer to "On the choice of step size in subgradient optimization"
"""
function get_stepsize!(state::DsState, δt_vec, δb_vec, L_obj::Float64)
    # theta_ub = ...
    (;r,r1, eps0, theta0, theta_c, beta_k) = state
    ar::Float64 = exp(-0.6933* ((r/r1)^3.26))
    if ar < eps0
        ar = eps0
        state.phaseI = false
    end
    ar = max(ar, eps0)
    theta_bar = ar*theta0 + (1-ar)*theta_c
    # theta_bar = theta0
    # div = sum(abs.(δt_vec)) + sum(abs.(δb_vec) )
    # div = sum(δt_vec .^ 2) + sum(δb_vec .^ 2)

    # Note: we divide the norm later in the get_new_dual function
    stepsize = 1/beta_k * (theta_bar - L_obj)
    # stepsize = 1/sqrt(beta_k)

    # Note: since we are smoothing the solution, D(λ) is not solved exactly, so it is possible that stepsize is negative. We may increase the obj upper bound with the hope that it is always positive
    # @assert (theta0 >= L_obj)
    if stepsize < 0
        # @warn "Got negative step size!" stepsize theta_bar L_obj
        stepsize = 1/sqrt(beta_k)
        # state.theta0 = L_obj
    end
    # stepsize = 1 / state.k 
    # @debug "Getting stepsize $stepsize" theta_bar L_obj ar
    return stepsize
end

"""
Get modified subgradient for better convergence.

Refer to the lecture notes by Boyd "Subgradient Method"
"""
function get_modified_grad!(state::DsState, δt_vec, δb_vec)
    if isempty(state.grad_b_m)
        state.grad_b_m = δb_vec
        state.grad_t_m = δt_vec
    else
        (;grad_b_m, grad_t_m) = state
        gamma = 0.1
        b =  - gamma * (sum(δt_vec.*grad_t_m) + sum(δb_vec.*grad_b_m))
        b /= (LA.norm(grad_b_m)^2 + LA.norm(grad_t_m)^2)
        b = max(0, b)
        state.grad_b_m = δb_vec + b * grad_b_m
        state.grad_t_m = δt_vec + b * grad_t_m
        # state.grad_b_m = δb_vec
        # state.grad_t_m = δt_vec 
    end
    return state.grad_t_m, state.grad_b_m 
end

function get_new_dual(dual::DsDual, stepsize::Float64, δt_vec, δb_vec)
    N = length(dual.λ) - 1
    dual1 = DsDual(N)
    b_div = sum(δb_vec.^2) + sum(δt_vec.^2)
    t_div = sum(δt_vec.^2) + sum(δb_vec.^2)

    for i in 1:N+1
        if dual.λ[i] > 0 || (δt_vec[i] > 0)
            # dual1.λ[i] = dual.λ[i] + stepsize / t_div * δt_vec[i]
        end
        dual1.λ[i] = dual.λ[i] + stepsize / t_div * δt_vec[i]

        if dual.μ[i] > 0 || (δb_vec[i] > 0)
            # dual1.μ[i] = dual.μ[i] + stepsize / b_div * δb_vec[i]
        end
        dual1.μ[i] = dual.μ[i] + stepsize / b_div * δb_vec[i]

        dual1.λ[i] = max(0.0, dual1.λ[i])
        dual1.μ[i] = max(0.0, dual1.μ[i])
    end
    return dual1
end

