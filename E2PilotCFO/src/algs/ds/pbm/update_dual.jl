
using MadNLP
using Ipopt

"""
Update the dual variable by solving a QP.
Let us use all the cutting planes for now...

return: state, the new dual variable
"""
function pbm_update_dual!(state::PbmState, primal, dual, δt_vec, δb_vec, L_obj::Float64, verbose::Int = 0)
    state.k += 1
    state.cons_vio = max(maximum(pos(δt_vec)),  maximum(pos(δb_vec)))
    (; uk, dual_center) = state

    lam_bar_vec = ds_dual2vec(dual_center)
    N = get_N(primal)

    if isinf(state.D_center)
        state.D_center = L_obj
    end

    push!(state.dual_vec, dual)
    push!(state.D_vec, L_obj)
    push!(state.g_vec, vcat(δt_vec, δb_vec))
    if L_obj >= state.kappa * state.delta + state.D_center
        state.D_center = L_obj
        state.dual_center = dual
    end

    if L_obj > state.theta_c
        # If there is an improvement 
        state.theta_c = L_obj 
        state.dual_c = deepcopy(dual)
        state.primal_c = deepcopy(primal)
    end

    model = Model(()->MadNLP.Optimizer(print_level=MadNLP.INFO,))
    # model = Model(OSQP.Optimizer)
    # model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, d)
    n_lam = 2*(N+1)
    @variable(model, lam_vec[1:n_lam] >= 1e-6)
    @objective(model, Max, d - 0.5 * uk *  sum( x^2 for x in (lam_vec - lam_bar_vec) ))
    for i in 1:length(state.g_vec)
        D = state.D_vec[i]
        dual = state.dual_vec[i]
        g = state.g_vec[i]
        lam_vec_i = ds_dual2vec(dual)
        @constraint(model, d <= D + sum( g.* (lam_vec .- lam_vec_i)), base_name="cons$i")
    end

    optimize!(model)
    lam_vec_new = value.(lam_vec)
    dual_new = ds_vec2dual(lam_vec_new)
    model_obj = objective_value(model)
    state.delta = model_obj  - state.D_center
    if verbose > 1
        @show model_obj state.D_center state.delta
    end
    
    return state, dual_new 
end