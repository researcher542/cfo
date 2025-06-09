
"""
Update the parameters for computing the stepsize. And update the dual variable.

return: the new dual variable.
"""
function update_state!(state::DsState, primal::DsPrimal, dual::DsDual, δt_vec, δb_vec, L_obj::Float64)

    state.k += 1

    state.cons_vio = max(maximum(pos(δt_vec)),  maximum(pos(δb_vec)))
    # update the step size
    if state.phaseI
        # @debug "Entering phase I"
        theta_c = state.theta_c
        cmp = isinf(theta_c) ? theta_c : theta_c + state.eps * abs(theta_c)
        if L_obj > cmp
            # If there is an improvement 
            state.theta_c = L_obj 
            state.dual_c = deepcopy(dual)
            state.primal_c = deepcopy(primal)
            state.v = 0
        else
            state.v += 1
            if state.v >= state.v_ub
                state.v = 0
                state.r += 1
            end
        end
    else
        # Phase II
        # @debug "Entering phase II"
        if L_obj > state.theta_c
            # If there is an improvement 
            state.theta_c = L_obj 
            state.dual_c = deepcopy(dual)
            state.primal_c = deepcopy(primal)
            state.v = 0
        end
        state.v += 1
        if state.v >= state.v_ub
            state.v = 0
            state.beta_k += 2
            # Reset to current best
            if state.beta_k < state.beta_ub
                dual = state.dual_c
            end
        end
    end

    stepsize = get_stepsize!(state, δt_vec, δb_vec, L_obj)
    state.stepsize = stepsize
    grad_t_m, grad_b_m = get_modified_grad!(state, δt_vec, δb_vec)
    dual = get_new_dual(dual, stepsize, grad_t_m, grad_b_m)
    return state, dual 
end