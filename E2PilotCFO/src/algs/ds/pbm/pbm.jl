
include("state.jl")
include("update_dual.jl")

"""
perform cfo_kpath and append the solution to it.
"""
# function append_kpath!(prob::AbsCfoProb, state0, option0)
# 
#     kp_option = KpOption(;K=16, restrict_N = true, fast_flag=false, thread_flag=true)
#     state_vec = cfo_kpath(prob, option0, kp_option)
#     kpath_state = state_vec[end]
# 
#     if kpath_state.obj_f < state0.obj_f
#         state0.primal_f = kpath_state.primal_f
#         state0.obj_f = kpath_state.obj_f
#         state0.cons_vio_f = kpath_state.cons_vio_f
#         state0.theta_f = kpath_state.theta_f
#     end
#     # @assert false
# 
#     return state0
# end

"""
The proximal bundle method that tries to solve the dual problem
"""
function cfo_bundle_method(prob::AbsCfoProb, option0)
    t_begin = time()
    option = copy(option0)
    (;maxiter, maxtime) = option
    state = PbmState()
    ds_preprocess!(prob, state, option)

    if option.early_break
        return [state]
    end
    N = option.N

    dual = DsDual(option.N)
    primal::DsPrimal = option.primal0
    state.dual_center = dual
    if isempty(primal)
        primal = DsPrimal(;sub_sol_vec = [DsSubsolution() for i in 1:N+1])
    end

    # if option.kpath_flag
    #     @timeit g_to "kpath" state = append_kpath!(prob, state, option)
    # end
    

    state_vec = PbmState[]
    ext_net = ds_construct_stage_ext_graph(prob, option)
    primal.tau_cs_vec = zeros(N+2)
    primal.beta_cs_vec = zeros(N+2)
    while state.k <= maxiter
        (option.verbose > 0) && println("iter=$(state.k) with obj=$(prob.objtype) ----------------")

        @timeit g_to "ds_solve_primal_ex" primal = ds_solve_primal_ext(prob, ext_net, primal, dual, option,)
        δt_vec, δb_vec, tau_cs_vec, beta_cs_vec, obj = get_subgradient(prob, primal, option)
        primal.tau_cs_vec = tau_cs_vec
        primal.beta_cs_vec = beta_cs_vec

        L_obj = get_lag_obj(obj, dual, δt_vec, δb_vec)
        if L_obj > (state.obj_f * (1+1e-2))
            @warn "D(λ) > f(x)!" L_obj state.obj_f prob.objtype 
            @show obj dual δt_vec δb_vec tau_cs_vec primal
            @exfiltrate
            @assert false "D(λ) > f(x)!"
        end

        # If termination condition is met
        if ( sum(δt_vec .^ 2) + sum(δb_vec .^ 2) < 1e-6 )
            break
        end
        if abs(state.delta) < 1e-6 * abs(state.obj_f)
            break
        end

        # update the debug information
        if option.verbose > 0
            @show (tau_cs_vec * option.t_mul) (beta_cs_vec * option.b_mul)
            @show δt_vec δb_vec
            @show dual_obj = (L_obj - obj)
        end

        state.obj = obj
        state.lag_obj = L_obj
        state.dual = dual
        state.primal = primal

        push_reopt_candidate!(state, primal)

        @timeit g_to "update dual" state, dual = pbm_update_dual!(state, primal, dual, δt_vec, δb_vec, L_obj, option.verbose)
        # @timeit g_to "update_primal" update_primal!(prob, state, primal, option)

        tt_ep = (time() - t_begin)
        if option.verbose > 0
            pbm_show_state(state) 
            @show tt_ep
            ds_print_debug(prob, primal, dual, option)
            println("\n")
        end

        if option.trace
            push!(state_vec, deepcopy(state))
        end

        if tt_ep > maxtime
            @warn "$tt_ep seconds excede maximum time " 
            if !isinf(state.obj_f)
                break
            end
        end
       
    end

    @timeit g_to "reopt" state = reopt_primal_vec!(prob, state, option)
    if option.trace
        push!(state_vec, deepcopy(state))
    end
    # @timeit g_to "update_primal" update_primal!(prob, state, primal, option)

    if option.trace
        return state_vec
    else
        return state
    end
end


function pbm_show_state(state)
    @show state.cons_vio 
    @show state.cons_vio_f
    @show state.lag_obj
    @show state.obj_f state.theta_c
    @show state.obj_f0
end