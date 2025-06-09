"""
Using the JuMP.jl to model 
"""

# import Juniper
import Gurobi

const GRB_ENV = Ref{Gurobi.Env}()
function init_gurobi()
    GRB_ENV[] = Gurobi.Env()
end

using PiecewiseLinearOpt
using Polynomials
using Infiltrator
using Dates
import SCIP


@with_kw mutable struct MIPOption
    t_mul::Float64 = 1/60.0
    b_mul::Float64 = 1/3.6e6
    verbose::Int = 0
    apx_type::Symbol = :ApxPiecewiseLinear
    time_limit::Float64 = 60.0 
    mem_limit::Float64 = 16.0 * 1024.0
end

"""
We are given some datapoints (x,y) and we want to find use a nonlinear function that approximates the data points.
"""
function apx_nlp_func(model, apx_type::Symbol, x_var, x_range, func; npoly::Int = 3)

    if length(x_range) < 2
        return func(x_range[1])
    elseif abs(x_range[end] - x_range[1]) < 1e-6
        return func(x_range[1])
    end


    if apx_type == :ApxPiecewiseLinear
        # # @info "Piecewise linear approximation"
        method =[:Logarithmic, :CC, :MC, :SOS2, :Incremental, 
              # default
            :DisaggLogarithmic, :ZigZag, :ZigZagInteger][1]
        return piecewiselinear(model, x_var, x_range, func; method=method) 
    elseif apx_type == :ApxPolynomial
        y_range = [func(x) for x in x_range]
        p = Polynomials.fit(x_range, y_range, npoly)
        val = sum(x_var^(i-1) * p.coeffs[i] for i in 1:length(p.coeffs))
        return val
    end
    @error("Unknown approximation type: $apx_type")
    return x_var 
end

function get_model_gurobi(option, region)
    model = Model(Gurobi.Optimizer)
    # set_optimizer_attribute(model, "NonConvex", 2) # We have bilinear terms.
    set_optimizer_attribute(model, "NumericFocus", 3) # To avoid numerical issues.
    set_optimizer_attribute(model, "TimeLimit", option.time_limit) # 
    set_optimizer_attribute(model, "SoftMemLimit", option.mem_limit / 1024) # 
    set_optimizer_attribute(model, "MemLimit", option.mem_limit / 1024) # 
    if option.verbose == 0
        set_optimizer_attribute(model, "LogToConsole", 0)
        # set_optimizer_attribute(model, "OutputFlag", 0)
    end

    log_dir = joinpath(ep.k_root_dir, "log", "$(today())", "grb")
    log_path = joinpath(log_dir, "grb-$(region)-$(now()).log")
    if !isdir(log_dir)
        try 
            mkdir(log_dir)
        catch e
            @warn "failed to create log dir $log_dir with error $e. But we will try to continue." 
        end
    end
    set_optimizer_attribute(model, "LogFile", log_path)

    # set_optimizer_attribute(model, "MIPFocus", 1) # Find more feasible solutions.
    # set_optimizer_attribute(model, "ScaleFlag", 2) # More aggressive scaling.
    # set_optimizer_attribute(model, "FeasibilityTol", 1e-9) # Tolerant more on the constraint violation, default is 1e-6. In the default value, Gurobi might have trouble in finding the feasible solution.
    nthreads = 8
    set_optimizer_attribute(model, "Threads", nthreads)  
    

    return model
end

function get_model_scip(option)
    @info "Getting SCIP solver model...."
    model = Model(SCIP.Optimizer)
    set_attribute(model, "display/verblevel", 3)
    set_attribute(model, "limits/gap", 0.01)
    set_attribute(model, "limits/time", option.time_limit)
    # set_attribute(model, "limits/time", 60.0)
    return model
end


cfo.cfo_mip_jump(prob, option::MIPOption) = cfo_mip_jump(prob, option)

"""
Given a initial primal solution, set the primal solution to the model.

Currently, we only consider the path
"""
function cfo_mip_jump_set_primal!(prob, model, primal, option)
    if isnothing(primal) || isempty(primal.sub_sol_vec)
        @debug "No primal solution provided, skip setting primal."
        return
    end
    (;net, ev, src, des) = prob

    (; t_mul, b_mul, verbose) = option

    selected_cs_vec = []
    for (isol, sol) in enumerate(primal.sub_sol_vec)
        if isol > 1
            nd = sol.path[1]
            @constraint(model, model[:y][nd] == 1)
            push!(selected_cs_vec, nd)
        end
        for i in 1:length(sol.path)-1
            u = sol.path[i]
            v = sol.path[i+1]
            @constraint(model, model[:x][u, v] == 1)
            @constraint(model, model[:t][u, v] == sol.t_vec[i] * t_mul)
        end 
    end

    cs_idx_vec = [cs.idx for cs in net.cs_vec if (cs.idx != prob.src && (cs.idx != prob.des))]
    for cs_idx in cs_idx_vec
        if !(cs_idx in selected_cs_vec)
            @constraint(model, model[:y][cs_idx] == 0)
            @constraint(model, model[:tc][cs_idx] == 0)
            @constraint(model, model[:Fhat][cs_idx] == 0)
        end
    end
    
end

"""
solve the cfo problem with MIP modelling and JuMP
"""
function cfo_mip_jump(prob, option = MIPOption(); primal0 = nothing)
    (;net, ev, src, des, T) = prob
    (; t_mul, b_mul, verbose) = option
    B = ev.cap
    B2 = 2*B
    n_node = ep.nv(net)
    cs_idx_vec = [cs.idx for cs in net.cs_vec if (cs.idx != prob.src && (cs.idx != prob.des))]
    # model = get_model_scip(option)
    model = get_model_gurobi(option, net.region)
    
    #### 
    way_vec = collect(keys(net.waydata))

   
    if verbose == 0
        set_silent(model)
    end
    @variable(model, y[1:n_node], Bin)
    
    @variable(model, x[way_vec], Bin)

    @variable(model, t[way_vec])

    @variable(model, 0 <= beta[1:n_node] <= B * b_mul)
    @variable(model, 0 <= beta_hat[1:n_node] <= B * b_mul) # The SoC after charging 

    @variable(model, 0 <= beta_bar[way_vec] <= B * b_mul) # The SoC after running on an edge

    @variable(model, 0 <= tau[1:n_node] <= T * t_mul)

    @variable(model, -B*b_mul <= c[way_vec] <= B * b_mul) # The energy consumption on a road segment

    @exfiltrate

    @variable(model, 0 <= Fhat[1:n_node])
    @variable(model, 0 <= tc[1:n_node] <= cfo.g_max_charge_time * t_mul)
    @variable(model, cfo.g_min_wait_time * t_mul <= tw[1:n_node]  <= cfo.g_max_wait_time * t_mul)

    @constraint(model, beta[src] == prob.β0 * b_mul)
    @constraint(model, tau[src] == 0.0)

    cfo_mip_jump_set_primal!(prob, model, primal0, option)

    for inode in 1:n_node
        outnei = outneighbors(net.g, inode)
        innei = inneighbors(net.g, inode)

        # Forbidden decisions that is not on the charging station.
        if !(inode in cs_idx_vec)
            @constraint(model, y[inode] == 0)
            @constraint(model, tc[inode] == 0)
            @constraint(model, Fhat[inode] == 0)
        end

        if (inode in cs_idx_vec)
            ## only charge on the selected path.
            @constraint(model, y[inode] <= sum(x[(inode, v)] for v in outnei))

            ## After charging, the SoC should be larger than the current SoC.
            @constraint(model, beta[inode] <= beta_hat[inode])

            ## Ensure that we charge only when we select it
            # @constraint(model, tc[inode] * (1-y[inode]) == 0) 
            @constraint(model, tc[inode] <= y[inode] * cfo.g_max_charge_time * t_mul)
            # for v in outnei
            #     ## only charge on the selected path.
            # end
        end

        # network flow 
        if inode == src
            flow_flag = 1
        elseif inode == des
            flow_flag = -1
        else
            flow_flag = 0
        end
        @constraint(model, sum(x[(inode, w)] for w in outnei) - sum(x[(u, inode)] for u in innei)  == flow_flag )

    end

    ## constraints for tau and beta
    for edge in edges(net.g)
        (u,v) = edge.src, edge.dst

        # travel time constraints 
        (min_t, max_t) = cfo.get_minmax_t(net, u, v, ev)
        ## Inverval constraint does not support conflict status
        @constraint(model, min_t * t_mul <= t[(u,v)] )
        @constraint(model, t[(u,v)] <= max_t * t_mul)

        # coupling constraints for tau
        if u in cs_idx_vec
            mid_var = tau[v] - tau[u] - t[(u,v)] - (tc[u] + tw[u]) * y[u]
            @constraint(model, -T * t_mul * (1-x[(u,v)]) <= mid_var)
            @constraint(model, mid_var <= t_mul * T*(1-x[(u,v)]))
            # @constraint(model,  -T * (1-x[(u,v)]) - T * (1 - x[u, v])  <=  <= T * (1 - x[u, v]))
        else
            mid_var = tau[v] - tau[u] - t[(u,v)] 
            @constraint(model, -T * t_mul * (1-x[(u,v)]) <= mid_var)
            @constraint(model, mid_var <= t_mul * T*(1-x[(u,v)]))
            # @constraint(model,  - T * (1 - x[u, v])  <= tau[v] - tau[u] - t[(u,v)] <= T * (1 - x[u, v]))
        end

        # coupling constraints for beta
        n_piece_ct = 10
        ct_func = t -> cfo.energy_cost_on_road(net, u, v, ev, t / t_mul) * b_mul
        min_t, max_t = cfo.get_minmax_t(net, u, v, ev)
        t_range = LinRange(min_t, max_t, n_piece_ct) * t_mul
        ct = apx_nlp_func(model, option.apx_type, t[(u,v)], t_range, ct_func)
        @constraint(model, c[(u,v)] == ct)

        n_piece_charge = 20
        if u in cs_idx_vec && (u != src) && (v != des) && (u != des)
            cf_inv_func = beta -> cfo.cf_inv(0.0, beta / b_mul, B) * t_mul
            cf_func = tc -> charge_function(0.0, tc / t_mul, B) * b_mul
            beta_range = collect(LinRange(0.0, B * b_mul, n_piece_charge))
            max_tc = cf_inv(0.0, B, B) * 1.2 * t_mul
            min_tc = g_min_charge_time * t_mul
            t_range = LinRange(min_tc, max_tc, n_piece_charge)
            tc0 = apx_nlp_func(model, option.apx_type, beta[u], beta_range, cf_inv_func)
            phi = apx_nlp_func(model, option.apx_type, tc0 + tc[u], t_range, cf_func)

            ## the final beta after charging.
            @constraint(model, phi == beta_hat[u])
            mid_var = beta[v] - beta_hat[u] + c[(u,v)]
            @constraint(model, - B2 * b_mul * (1-x[(u,v)]) <= mid_var)
            @constraint(model,  mid_var <= B2 * b_mul * (1-x[(u,v)]))
        else
            @constraint(model, beta_bar[(u,v)] <= (beta[u] - c[(u,v)]) )
            mid_var = beta[v] - beta_bar[(u,v)]
            @constraint(model, - B2 * (1-x[(u,v)]) * b_mul <= mid_var)
            @constraint(model,  mid_var <= B2 * (1-x[(u,v)]) * b_mul)
        end

    end

    ## constraints for objective
    
    for inode in cs_idx_vec
        t_range = collect(0:3600:(T+3600.0)) * t_mul
        ta = get_price_ta(prob, inode, prob.predict_mode)
        pi_func = tau -> get_price(ta, prob.start_time, tau / t_mul)[1]

        pi_v = apx_nlp_func(model, option.apx_type, tau[inode] + tw[inode], t_range, pi_func)
        Fv = pi_v * (beta_hat[inode] - beta[inode])
        @constraint(model, Fhat[inode] >= Fv - b_mul * B * 2.0 * (1-y[inode]) )
    end
    init_ci = cfo.get_init_battery_ci(prob, prob.predict_mode)

    F_end = init_ci * (prob.β0 * b_mul - beta[des])

    @show init_ci
    ## 
    @objective(model, Min, sum(Fhat) + F_end )

    

    @info "Start solving the model with JuMP.jl"
    JuMP.optimize!(model)

    status = termination_status(model)
    if verbose > 0
        solution_summary(model)
    end

    primal = model2primal(prob, model, option)

    ## convert the model to a DsPrimal solution

    return primal, model
end

"""
Convert the solved model to a DsPrimal solution

return: DsPrimal
"""
function model2primal(prob, model, option)
    status = termination_status(model)
    @info "Converting the model to primal solution with status $(status)"
    (;net, src, des) = prob
    (; t_mul, b_mul) = option
    if status == INFEASIBLE
    # iis_model, vec = cfo.debug_infeasible_solution(model)
        @warn "model is infeasible."
        if false
            vec = debug_infeasible_solution(model)
            for cons in vec
                println(cons)
                println()
            end
        end
        return DsPrimal()
    elseif primal_status(model) == NO_SOLUTION
        @warn "No solution found for the model."
        return DsPrimal()
    end
    x_mat = value.(model[:x]) 
    y_vec = value.(model[:y])
    t_mat = value.(model[:t]) / t_mul
    tc_vec = value.(model[:tc]) / t_mul
    tw_vec = value.(model[:tw]) / t_mul
    F_vec = value.(model[:Fhat])
    beta_hat_vec = value.(model[:beta_hat]) 
    beta_vec = value.(model[:beta]) 
    tau_vec = value.(model[:tau]) 

    F_vec1 = filter(x->x!=0, F_vec)
    @show F_vec1

    sub_sol_vec = DsSubsolution[]
    sub_path = Int[src]
    sub_t_vec = Float64[]
    tc = tc_vec[src]
    tw = tw_vec[src]
    if y_vec[src] == 0
        tc = tw = 0.0
    end
    if tc == 0.0
        tw = 0.0
    end
    cnt = 0
    @show length(y_vec)
    cur_node = -1
    while length(sub_sol_vec) < 1000 && length(sub_path) <= length(y_vec) && cnt <= length(y_vec)
        # @info "sub_path: $(sub_path), sub_t_vec: $(sub_t_vec), cnt: $cnt, cur_node: $(sub_path[end])"
        cnt += 1
        if cur_node == sub_path[end]
            @warn "The current node is the same as the last node, break the loop. Something wrong might happen."
            break
        end
        cur_node = sub_path[end]
        if y_vec[cur_node] == 1 && cur_node != des
            nd = sub_path[1]
            tc = tc_vec[nd]
            tw = tw_vec[nd]
            if tc == 0.0 && nd == src
                tw = 0.0
                tw_vec[nd] = 0.0
            end
            if option.verbose > 0
                @show F_vec[nd]
                @show beta_vec[nd] beta_hat_vec[nd] tc/60
                @show tw / 60.0
                Δβ = beta_hat_vec[nd] - beta_vec[nd]
                @show Δβ 
                @show tau_vec[nd]
                println("*"^80)
            end
            ## If it is the selected charging station.
            sub_sol = DsSubsolution(path=sub_path, t_vec=sub_t_vec, tc=tc, tw=tw)
            if cur_node != src
                push!(sub_sol_vec, sub_sol)
            end
            sub_path = Int[cur_node]
            sub_t_vec = Float64[]
            ## The charging and waiting time for the next stage.
            # tc = tc_vec[cur_node]
            # tw = tw_vec[cur_node]
        end
        if cur_node == des
            if !isempty(sub_t_vec)
                nd = sub_path[1]
                path_beta_vec = value.(model[:beta][sub_path])
                # path_beta_bar_vec = value.(model[:beta_bar][sub_path])
                path_tau_vec = value.(model[:tau][sub_path])
                path_ct_vec = value.([model[:c][(sub_path[i], sub_path[i+1])] for i in 1:length(sub_path)-1])
                tc = tc_vec[nd]
                path_t_vec = value.([model[:t][(sub_path[i], sub_path[i+1])] for i in 1:length(sub_path)-1])
                tw = tw_vec[nd]
                if tc == 0.0 && nd == src
                    tw = 0.0
                    tw_vec[nd] = 0.0
                end
                if option.verbose > 0
                    @show F_vec[nd]
                    @show beta_vec[nd] beta_hat_vec[nd] tc/60
                    Δβ = beta_hat_vec[nd] - beta_vec[nd]
                    @show tw / 60.0
                    @show Δβ 
                    @show tau_vec[nd]
                    @show tau_vec[cur_node]
                    # @show path_beta_vec[1:15]
                    # @show path_beta_bar_vec[1:15]
                    #@show path_tau_vec[1:15]
                    # @show path_t_vec[1:15]
                    # @show path_ct_vec[1:15]
                    @show length(path_beta_vec)
                    println("*"^80)
                end
                sub_sol = DsSubsolution(path=sub_path, t_vec=sub_t_vec, tc=tc, tw=tw)
                push!(sub_sol_vec, sub_sol)
            end
            break
        end
        for v in outneighbors(net.g, cur_node)
            # @info "Checking node $cur_node to $v, x_mat[cur_node, v] = $(x_mat[cur_node, v])"
            if abs(x_mat[(cur_node, v)] - 1) < 1e-4
                push!(sub_path, v)
                push!(sub_t_vec, t_mat[(cur_node, v)])
                break
            end
        end
    end
    if isempty(sub_sol_vec)
        @warn "No subsolution found, return empty DsPrimal."
        @exfiltrate
        throw(ErrorException("No subsolution found, return empty DsPrimal."))
    end
    primal = DsPrimal(sub_sol_vec=sub_sol_vec)
    return primal
end


function debug_infeasible_solution(model::Model)
    status = termination_status(model)

    list_of_conflicting_constraints = ConstraintRef[]
    if status != INFEASIBLE
        @warn "model is feasible with status $status."
        return list_of_conflicting_constraints
    end

    compute_conflict!(model)
    for (F, S) in list_of_constraint_types(model)
        if S in [MOI.ZeroOne] || F in [VariableRef]
            ## Interval does not support conflict status
            # continue
        end
        for cons in all_constraints(model, F, S)
            if get_attribute(cons, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
              push!(list_of_conflicting_constraints, cons)
              # @show cons F S
              # @show cons
            end
        end
    end
    return list_of_conflicting_constraints
    
    # iis_model, _ = copy_conflict(model)
    # for (F, S) in list_of_constraint_types(iis_model)
    #     for con in all_constraints(iis_model, F, S)
    #         push!(list_of_conflicting_constraints, con)
    #     end
    # end
    # return iis_model, list_of_conflicting_constraints
end