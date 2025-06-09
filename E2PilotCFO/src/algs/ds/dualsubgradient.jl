"""
The dual subgradient method. 
we use ds_ as its prefix for related functions
"""

include("primal_var.jl")
include("option.jl")
include("dual_var.jl")


function cfo_dual_subgradient(prob0::AbsCfoProb, option0::DsOption)
    t_begin = time()
    option = deepcopy(option0)
    prob = get_cfo_prob_fix(prob0)
    (;maxiter, maxtime) = option
    state = DsState()
    @timeit g_to "ds_preprocess" ds_preprocess!(prob, state, option)
    if option.early_break
        return [state]
    end
    N = option.N
    if option.verbose > 0
        @debug "" option.N
    end

    dual = DsDual(option.N)

    primal = DsPrimal(;sub_sol_vec = [DsSubsolution() for i in 1:N+1])

    state_vec = DsState[]

    ext_net = ds_construct_stage_ext_graph(prob, option)
    primal.tau_cs_vec = zeros(N+2)
    primal.beta_cs_vec = zeros(N+2)
    while state.k <= maxiter
        (option.verbose > 0) && println("iter=$(state.k) with obj=$(prob.objtype) ----------------")

        @timeit g_to "ds_solve_primal" primal = ds_solve_primal_ext(prob, ext_net, primal, dual, option,)
        # @timeit g_to "ds_solve_primal" primal = ds_solve_primal(prob, primal, dual, option,); # @exfiltrate; @assert false
        δt_vec, δb_vec, tau_cs_vec, beta_cs_vec, obj = get_subgradient(prob, primal, option)
        primal.tau_cs_vec = tau_cs_vec
        primal.beta_cs_vec = beta_cs_vec
        L_obj = get_lag_obj(obj, dual, δt_vec, δb_vec)

        if L_obj > (state.obj_f * (1+1e-2)) && option.verbose > 0
            @warn "D(λ) > f(x)!" L_obj state.obj_f prob maxlog=3
            # @assert false "D(λ) > f(x)!"
        end

        # If termination condition is met
        if ( sum(δt_vec .^ 2) + sum(δb_vec .^ 2) < 1e-6 )
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
        state, dual = update_state!(state, primal, dual, δt_vec, δb_vec, L_obj)
        # @timeit g_to "update_primal" update_primal!(prob, state, primal, option)
        tt_ep = (time() - t_begin)
        if option.verbose > 0
            ds_show_state(state) 
            @show tt_ep
            ds_print_debug(prob, primal, dual, option)
            println("\n")
        end

        if option.trace
            push!(state_vec, deepcopy(state))
        end

        if tt_ep > maxtime
            @warn "$tt_ep seconds exceede maximum time iter=$(state.k). break." 
            if !isinf(state.obj_f)
                break
            end
        end
       
    end

    @info "reoptimizing the primal paths..."
    @timeit g_to "reopt_primal_vec" state = reopt_primal_vec!(prob, state, option)
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

function ds_show_state(state)
    @show state.stepsize 
    @show state.cons_vio 
    @show state.cons_vio_f
    @show state.lag_obj
    @show state.obj_f state.theta_c
    @show state.obj_f0
end

function get_lag_obj(prob::AbsCfoProb, primal::DsPrimal, dual::DsDual, option)
    δt_vec, δb_vec, tau_cs_vec, beta_cs_vec, obj = get_subgradient(prob, primal, option)
    return get_lag_obj(obj, dual, δt_vec, δb_vec) 
end

function get_lag_obj(obj::Float64, dual::DsDual, δt_vec, δb_vec, debug::Bool=false) 
    at = sum(dual.λ .* δt_vec)
    ab = sum(dual.μ .* δb_vec)
    if debug
        @show obj at ab dual
        @show δb_vec δt_vec
    end
    return obj +  at + ab
end

function ds_print_debug(prob::AbsCfoProb, primal::DsPrimal, dual::DsDual, option::DsOption)
    (;t_mul, b_mul) = option
    sol_vec = primal.sub_sol_vec
    sol_end = primal.sub_sol_vec[end]
    des = sol_end.path[end]
    cs_path = vcat([sol.path[1] for sol in primal.sub_sol_vec], des)
    τ_des = sol_end.τ + sol_end.tw + sol_end.tc + sum(sol_end.t_vec)
    tau_vec_cs = vcat([sol.τ for sol in primal.sub_sol_vec], τ_des) * t_mul
    beta_vec_cs= [sol.β for sol in primal.sub_sol_vec] * b_mul
    tc_vec = [sol.tc for sol in sol_vec]
    tw_vec = [sol.tw for sol in sol_vec]
    (;obj) = ds_simulate(prob, primal; check_flag=false, predict_mode=prob.predict_mode)
    @show obj 
    # @show beta_vec tau_vec
    # @show beta_vec_cs tau_vec_cs
    @show dual.λ dual.μ
    @show tc_vec/60.0 tw_vec/60.0
end


include("state.jl")
include("solveprimal/solveprimal.jl")
include("stepsize.jl")
include("sim.jl")
include("get_sub_grad.jl")
include("update_state.jl")
include("reopt/reopt.jl")
include("sigma/sigma.jl")
include("preprocess.jl")
include("fast.jl")

include("pbm/pbm.jl")
include("kpath/kpath.jl")