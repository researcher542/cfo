

"""
construct a dummy primal from the cs_path
"""
function ds_cs_path2primal(cs_path::Vector{Int})
    sub_sol_vec = [
        DsSubsolution(path=[cs_path[i], cs_path[i+1]])
        for i in 1:length(cs_path)-1
    ] 
    push!(sub_sol_vec[end].path, cs_path[end])
    primal = DsPrimal(;sub_sol_vec=sub_sol_vec)
    return primal
end

function ds_make_cs_path_feasible(prob::AbsCfoProb, primal::DsPrimal)
    if ds_check_cs_path_feasible(prob, primal)
        return primal
    end

    cs_path = get_cs_path_unique(primal)
    # @show cs_path
    cs_path_new = Int[]
    for i in 1:length(cs_path)-1
        u = cs_path[i]
        v = cs_path[i+1]
        if u == v
            continue
        end
        push!(cs_path_new, u)
        if u == v
            continue
        elseif !has_edge(prob.net.cs_net, u, v)
            function getdist0(net, u, v, lam)
                cs_edge = get_edge(net, u, v)
                if isnothing(cs_edge)
                    return Inf
                else
                    return cs_edge.min_t
                end
            end
            inner_cs_path = shortest_path(Astar(), prob.net.cs_net, u, v; getdist = getdist0 )
            cs_path_new = vcat(cs_path_new, inner_cs_path[2:end-1])
            # @debug "" u v cs_path_new inner_cs_path
        end
    end
    push!(cs_path_new, cs_path[end])

    primal_new = ds_cs_path2primal(cs_path_new) 
    return  primal_new
end

"""
Check if the cs_path of a primal is feasible with the help of cs_net.
"""
function ds_check_cs_path_feasible(prob::AbsCfoProb, primal::DsPrimal)
    cs_path = get_cs_path(primal)
    for i in 1:length(cs_path)-1
        u = cs_path[i]
        v = cs_path[i+1]
        if u == v
            continue
        elseif !has_edge(prob.net.cs_net, u, v)
            return false          
        end
    end
    return true
end

function push_reopt_candidate!(state, primal::DsPrimal)
    h = hash_path(primal)
    if h in state.path_hash_set
        return
    end
    push!(state.reopt_candidate_vec, deepcopy(primal))
    push!(state.path_hash_set, h)
    return state 
end

function reopt_primal_vec!(prob::AbsCfoProb, state, option::DsOption)
    function func(primal)
        return ds_reopt(prob, primal, option)
    end

    primal_vec = state.reopt_candidate_vec
    if option.kpath_flag
        kp_option = KpOption(;K=16, restrict_N = true, fast_flag=false, thread_flag=true)
        kpath_primal_vec = kpath_get_primal_vec(prob, option, kp_option)
        primal_vec = vcat(primal_vec, kpath_primal_vec)
    end
    if option.verbose > 0
        @info "reoptimizing the primal candidates num=$(length(primal_vec))"
        # for primal in primal_vec
        #     @debug "reoptimizing $(get_cs_path(primal)), length=$(length(primal.sub_sol_vec)) N=$(get_N(primal))"
        # end
    end


    if isempty(primal_vec)
        return state
    end
    GC.gc(true)

    res_vec = Any[nothing for i in 1:length(primal_vec)]
    Threads.@threads :greedy for i in eachindex(primal_vec)
    # @batch for i in eachindex(primal_vec)
        primal = primal_vec[i]
        if option.verbose > 0
            @debug "reoptimizing $(get_cs_path(primal)), length=$(length(primal.sub_sol_vec)) N=$(get_N(primal))"
        end
        res_vec[i] = ds_reopt(prob, primal, option)
    end
    
    # if option.verbose > 0
    #      res_vec = ep.e2map(func, primal_vec, true, true; prefix = "reopt_primal_vec: ")
    # else
    #     res_vec = ep.e2map(func, primal_vec)
    # end
    res_vec = sort!(res_vec, by=x->x[1])
    tol = 1e-1
    if option.verbose > 0
        @info "reoptimizing the primal candidates failed num=$(sum(isinf.([res[1] for res in res_vec])))"
        for (solver_obj1, primal1) in res_vec
            @debug "reoptimized $(get_cs_path(primal1)), length=$(length(primal1.sub_sol_vec)), solver_obj=$solver_obj1"
        end
    end
    for ires in 1:length(res_vec)
        res = res_vec[ires]
        solver_obj, primal1 = res
        if option.verbose > 0
            @debug "Got result with solver_obj=$solver_obj, length=$(length(primal1.sub_sol_vec))"
        end
        if isinf(solver_obj)
            if option.verbose > 0
                @debug "re-optimize failed with solver_obj=$solver_obj"
            end
            continue
        end

        debug_flag = false
        δt_vec, δb_vec, tau_cs_vec, beta_cs_vec, obj = get_subgradient(prob, primal1, option, debug_flag)
        cons_vio = max(maximum(pos(δt_vec)),  maximum(pos(δb_vec)))
        # @warn ires, cons_vio, obj, solver_obj
        if cons_vio > tol 
            ## It can happen that the objective is smaller and it is infeasible.
            # But we have map it back
            # @warn "reopt_primal_vec failed with cons_vio=$cons_vio solver_obj=$solver_obj obj=$obj δt_vec=$(δt_vec) δb_vec=$(δb_vec)"
            continue
        end
        if abs(solver_obj - obj) > abs(obj) * 2e-2 && !isinf(solver_obj)
            @warn "objective inconsistence in reopt_primal" solver_obj obj
        end
        if (cons_vio < tol) 
            if (state.obj_f < obj) && (state.cons_vio_f < tol)
                @debug "objective is larger than the stored value..." obj solver_obj state.obj_f
                continue
            end
            if option.verbose > 0
                @debug "re-optimize successfully with obj=$obj, solver objective=$solver_obj, previous=$(state.obj_f), cons_vio=$(cons_vio)."
            end
            state.primal_f = copy(primal1)
            state.theta_f = obj
            state.cons_vio_f = cons_vio
            state.obj_f = obj
            return state
        end
    end
    return state
end

"""
update and store the primal solution which is feasible.
"""
function update_primal!(prob::AbsCfoProb, state, primal::DsPrimal, option::DsOption)
    # (;alpha) =  state

    # check have we optimized this path
    h = hash_path(primal)
    if h in state.path_hash_set
        return
    end
    @debug "reoptimizing primal solution." get_cs_path(primal)
    push!(state.path_hash_set, h)

    solver_obj, primal_tmp = ds_reopt(prob, primal, option)

    if !isvalid(primal_tmp)
        @debug "re-optimize failed with solver_obj=$solver_obj, not valid primal."
        return
    end
    # if isinf(solver_obj)
    #     @debug "re-optimize failed with solver_obj=$solver_obj"
    #     return
    # end

    debug_flag = false
    δt_vec, δb_vec, tau_cs_vec, beta_cs_vec, obj = get_subgradient(prob, primal_tmp, option, debug_flag)
    cons_vio = max(maximum(pos(δt_vec)),  maximum(pos(δb_vec)))

    cons_tol = 1e-1


    # if !(solver_obj ≈ obj )
    #     @warn "solver obj and obj from get_subgradient not consistent." obj solver_obj
    # end

    lock(state.state_lock) do
        if (state.cons_vio_f > cons_tol && cons_vio < state.cons_vio_f)
            state.primal_f = copy(primal_tmp)
            state.theta_f = solver_obj
            @debug "re-optimize with less cons_vio, with obj=$obj, solver objective=$solver_obj, previous=$(state.obj_f), cons_vio=$(cons_vio), previous=$(state.cons_vio_f)."
            # @show δt_vec δb_vec
            state.obj_f = obj
            state.cons_vio_f = cons_vio

        elseif (obj < state.obj_f && cons_vio < cons_tol) 
            state.primal_f = copy(primal_tmp)
            state.theta_f = solver_obj
            state.cons_vio_f = cons_vio
            @debug "re-optimize successfully with obj=$obj, solver objective=$solver_obj, previous=$(state.obj_f), cons_vio=$(cons_vio)."
            # @show δt_vec δb_vec
            # @show tau_cs_vec beta_cs_vec
            state.obj_f = obj
        else
            @debug "re-optimize failed with cons_vio=$cons_vio or larger obj=$obj, previous=$(state.obj_f)" maximum(pos(δt_vec)) maximum(pos(δb_vec)) 
            reopt_δt_vec = δt_vec
            reopt_δb_vec = δb_vec
            # @show reopt_δt_vec reopt_δb_vec
        end
    end

end



