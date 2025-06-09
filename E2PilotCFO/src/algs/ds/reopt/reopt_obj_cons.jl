

function ds_reopt_cs_net_objective(t_vec_all, prob::AbsCfoProb, primal::DsPrimal, option::DsOption, debug::Bool = false)
    penalty = 1e10
    obj, b_cons_vec, t_cons = ds_reopt_cs_net_obj_cons(t_vec_all, prob, primal, option, debug)
    obj1 = obj + sum(pos.(b_cons_vec)) * penalty + pos(t_cons) * penalty
    return obj1
end

"""
return the objective and constraints all at once.
return (obj, b_cons_vec, t_cons)
"""
function ds_reopt_cs_net_obj_cons(t_vec_all, prob::AbsCfoProb, primal::DsPrimal, option::DsOption, debug::Bool = false)
    (;net, T, β0, ev) = prob
    (;cs_net) = net
    (;β_lb, t_mul, b_mul, N) = option
    N = get_N(primal)
    cs_path = get_cs_path(primal)
    t_vec, tc_vec, tw_vec = ds_reopt_cs_net_decompose(t_vec_all, primal) 

    #######################
    tau_vec = zeros(N+2)
    beta_vec = zeros(N+2)
    beta_vec[1] = β0
    obj = 0.0
    b_cons_vec = zeros(N+2)
    for istage in 2:N+2
        # tc = tc_vec[i]
        # tw = tw_vec[i]
        cs_idx = cs_path[istage-1]
        if cs_path[istage-1] == cs_path[istage]
            e_cost = 0.0 
            charged = 0.0
        else
            edge = get_edge(cs_net, cs_path[istage-1], cs_path[istage])
            e_cost = edge(t_vec[istage-1], true)
            # @timeit g_to "paso" _, _, e_cost = paso(net, cs_path[istage-1], cs_path[istage], t_vec[istage-1]; veh=ev, output_step = false, debug_msg=false, visitonce=false, breakearly=false)

            charged = charge_function(beta_vec[istage-1], tc_vec[istage-1], ev.cap) - beta_vec[istage-1]
        end
        beta_vec[istage] = beta_vec[istage-1] + charged - e_cost
        tau_vec[istage] = tau_vec[istage-1] + tc_vec[istage-1] + tw_vec[istage-1] + t_vec[istage-1]

        tau_vec[istage-1] =  clamp(tau_vec[istage-1], 0.0, Inf)
        tw_vec[istage-1] = clamp(tw_vec[istage-1], 0.0, Inf)
        obj1 = cfo_objective(prob, cs_idx, 
        beta_vec[istage-1], tc_vec[istage-1], 
        tau_vec[istage-1] + tw_vec[istage-1], 
        prob.predict_mode
        )

        # b_penalty += b_mul * pos(β_lb - beta_vec[istage])
        # b_cons = b_mul * pos(β_lb - beta_vec[istage])
        b_cons = b_mul * (β_lb - beta_vec[istage])
        # push!(b_cons_vec, b_cons)
        b_cons_vec[istage-1] = b_cons
        
        obj += obj1 * b_mul
        if debug
            @debug "At i=$(istage-1), tau=$(tau_vec[istage-1]) tr=$(t_vec[istage-1]) tw=$(tw_vec[istage-1]) tc=$(tc_vec[istage-1]) β0=$(beta_vec[istage-1]/3.6e6) obj=$(obj1/3.6e6) charged=$(charged/3.6e6) cs_idx=$cs_idx, e_cost=$(e_cost/3.6e6)"
        end
        if prob.objtype == ObjTime()
            obj += tw_vec[istage-1] + t_vec[istage-1]
        elseif prob.objtype == ObjEnergy() || prob.objtype == ObjCarbon()
            # add a small penalty for tw so that it does not wait too much with no benefit.
            obj += 1e-10 * (tw_vec[istage-1] + tc_vec[istage-1])
        end
        # @debug "" i beta_tmp charged
        # @show beta_vec
    end
    if prob.objtype == ObjTime()
        obj = tau_vec[end]
    end
    # b_penalty += b_mul * pos(β_lb - beta_vec[end])
    b_cons = b_mul * (β_lb - beta_vec[end])
    # push!(b_cons_vec, b_cons)
    b_cons_vec[end] = b_cons
    t_cons = t_mul * (tau_vec[end] - T) 
    # @show obj b_penalty t_penalty
    # obj_dual = obj + penalty * b_penalty + penalty * t_penalty
    # @show beta_vec*b_mul

    return obj, b_cons_vec, t_cons
end

 ##########################
    # Get the cost with paso and simulation.
    # primal_new = ds_reopt_cs_net_tvec2primal(t_vec_all, primal, prob, option)
    # obj, beta_vec, tau_vec = ds_simulate(prob, primal_new; check_flag = false)

    # a1 = obj
    # neg(v) = [x < 0 ?  x : 0.0 for x in v]
    # # sol_beta_vec = [sol.β for sol in primal_new.sub_sol_vec] 
    # beta_err_vec = [pos(β_lb - sol.β) for sol in primal_new.sub_sol_vec] * b_mul
    # a2 = penalty * sum(beta_err_vec)
    # a3 = penalty * pos(β_lb*b_mul - beta_vec[end]) 
    # a4 = penalty * sum([x < 0.0 ? -x : 0.0 for x in beta_vec])
    # a5 = 10 * penalty * pos(tau_vec[end] - T*t_mul)
    # # @show a1 a2 a3 a4 a5
    # obj_dual = a1 + a2 + a3 + a4 + a5
########################

