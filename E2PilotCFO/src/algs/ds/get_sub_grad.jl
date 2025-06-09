
"""
return: (δt_vec, δb_vec, tau_cs_vec, beta_cs_vec, obj)
Here obj is the objective without Lagrangian
"""
function get_subgradient(prob::AbsCfoProb, primal::DsPrimal, option::DsOption, debug::Bool = false)
    (;β_lb, b_mul, t_mul) = option
    N = get_N(primal)
    (;ev, net, β0, T) = prob
    B = ev.cap
    δt_vec = zeros(N+1)
    δb_vec = zeros(N+1)
    obj::Float64 = 0.0
    # The vectors that is the true information from "simulate" the solution.
    # This includes the src and des information.
    beta_cs_vec = zeros(N+2)
    beta_cs_vec[1] = β0
    tau_cs_vec = zeros(N+2)
    obj1::Float64 = 0.0
    for istage in 1:N+1
        sub_sol = primal.sub_sol_vec[istage]
        (;path, t_vec, tc, tw, τ, β, ics) = sub_sol
        next_tau = (istage == N+1) ? T : primal.sub_sol_vec[istage+1].τ


        δt_vec[istage]= sum(t_vec) + tw + tc - (next_tau - τ)
        tau_cs_vec[istage+1] = tau_cs_vec[istage] + sum(t_vec) + tw + tc

        # update for the soc constraint.
        Δβi = charge_function(β, tc, B) - β
        
        cost_vec = [energy_cost_on_road(net, path[i], path[i+1], ev, t_vec[i]) for i in 1:length(path)-1]

        beta_cs_vec[istage+1] = beta_cs_vec[istage] + Δβi - sum(cost_vec)
        if istage != 1
            # Note that here the τ should use the scheduled one instead of simulate one. Otherwise D(λ) is not actuallly the optimal
            obj1 = cfo_objective(prob, path[1], β, tc, τ + tw, prob.predict_mode)
            obj += obj1
        end
        
        next_β = (istage == N+1) ? β_lb : primal.sub_sol_vec[istage+1].β
        δb_vec[istage] = sum(cost_vec) - (β+Δβi) + next_β
        if debug 
            @debug "At i=$(istage), tau=$(tau_cs_vec[istage]) tr=$(sum(t_vec)) tw=$(tw), tc=$tc, β0=$(β/3.6e6) obj=$(obj1/3.6e6) charged=$(Δβi/3.6e6), e_cost=$(sum(cost_vec)/3.6e6)"
        end
        # @show i sum(cost_vec)/3.6e6 β/3.6e6 Δβi/3.6e6 next_β/3.6e6 ics δb_vec[i]/3.6e6
    end
    if prob.objtype == ObjTime()
        obj = tau_cs_vec[end] / b_mul * t_mul
    else
        used_init_beta = β0 - beta_cs_vec[end]
        if prob.objtype == ObjEnergy()
            init_ci = 1.0
        else
            # ta = get_price_ta(prob, prob.src, prob.predict_mode)
            # init_ci = get_price(ta, prob.start_time, 0.0)[1]
            init_ci = get_init_battery_ci(prob, prob.predict_mode)
        end
        ## accounts for the initial battery level
        obj += used_init_beta * init_ci
    end



    res = δt_vec*t_mul, δb_vec*b_mul, tau_cs_vec, beta_cs_vec, obj*b_mul
    return res
end

