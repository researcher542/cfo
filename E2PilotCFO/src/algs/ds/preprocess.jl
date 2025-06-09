"""
use fastest path on the cs_net to get some ideas of the problem instance.
"""
function ds_preprocess!(prob::AbsCfoProb, state, option::DsOption)
    (;net, src, des, ev, β0, T) = prob
    (;b_mul, t_mul) = option
    fast_step = fastest_path(net, src, des; veh=ev, ) 
    # slow_step = fastest_path(net, src, des; veh=ev, slow_flag = true) 
    # @show fast_step slow_step

    ds_update_fast_vec!(prob, option)

    
    t2 = Threads.@spawn begin
        primal2, flag2 = cfo_practice(prob)
        if !flag2
            # @timeit g_to "reopt_practice"  
            update_primal!(prob, state, primal2, option)
        end
    end

    t3 = Threads.@spawn begin
        if !isempty(option.primal0)
            ## If primal0 is non-empty, append it to the results.
            @debug "primal0 is non-empty, append it the the primal candidates."
            update_primal!(prob, state, option.primal0, option)
        end
    end

    t1 = Threads.@spawn begin
        @debug "Getting feasible path with cs_net..."
        primal1 = ds_cs_net_get_feasible_path(prob, option, false)
        update_primal!(prob, state, primal1, option)
    end


    @timeit g_to "reopt_feasible_path" wait(t1)
    @timeit g_to "reopt_practice" wait(t2)
    @timeit g_to "reopt_primal0" wait(t3)

    # @timeit g_to "get_practice" 
    # ## If it is feasible
    
    # @timeit g_to "get_feasible_path" 
    # @show primal1

    # @timeit g_to "reopt_feasible_path" 

    ## use another method (PRACTICE) to get a feasible path.
    
    @debug "cs path on cs_net with N=$(option.N)"

   
    if option.N == -1
        option.N = get_N(state.primal_f) + 2
    end
    state.obj_f0 = state.obj_f

    # state.theta0 = option.N * prob.ev.cap * option.b_mul
    # fast_cost = fast_step.summary.cost
    # slow_cost = slow_step.summary.cost
    # charged_energy = max(slow_cost-β0, 0.0)
    ## Note: fastest path need not to be the one with most of the energy, because we have to go to other charging stations and consume more energy...
    # state.theta0 = charged_energy * b_mul 
    # if prob.objtype == ObjTime()
    #     state.theta0 = prob.T * 3.0
    # else
    #     state.theta0 = option.N * ev.cap * b_mul
    # end

    ## A feasible objective is simply an upper bound of D(λ)
    state.theta0 = state.obj_f
    if T <= fast_step.summary.duration 
        @warn "T=$T is less than fastest path T=$(fast_step.summary.duration). One might need to modify the problem. We should also note the difference between fast and fast-ds. fast-ds is the true fastest path but requires more time, we do not include it here."
        # sleep(10)
        # error("")
    end
    
end

"""
update the minimum time from a src to each node
"""
function ds_update_fast_vec!(prob::AbsCfoProb, option::DsOption)
    (;net, ev, src, des) = prob
    I = [e.src for e in edges(net.g)]
    J = [e.dst for e in edges(net.g)]
    t_forward = [minmax_t(net, e.src, e.dst, ev)[1] for e in edges(net.g)]
    t_rev = t_forward
    if prob.net.continent == :us
        ## For the US map, we have that the forward and reverse time are different.
        t_rev = [minmax_t(net, e.dst, e.src, ev)[1] for e in edges(net.g)]
    end
    option.fast_t_vec = shortest_path_all(net, src, sparse(I, J, t_forward))
    option.fast_t_rev_vec = shortest_path_all(net, des, sparse(I, J, t_rev))
    return option
end


"""
update tau and beta of a primal variable based on the path and t_vec and tw and tc.
"""
function ds_update_tau_beta!(prob::AbsCfoProb, primal::DsPrimal, option::DsOption)
    (;t_mul, b_mul) = option
    (;obj, beta_vec, tau_vec) = ds_simulate(prob, primal; check_flag = false, predict_mode = prob.predict_mode)
    idx = 1
    for (isol,sol) in enumerate(primal.sub_sol_vec)
        sol.τ = tau_vec[idx] / t_mul
        sol.β = beta_vec[idx] / b_mul
        # sol.τ = min(sol.τ, prob.T)
        # sol.β = min(sol.β, option.β_lb)
        idx += length(sol.t_vec)
    end
    return primal
end