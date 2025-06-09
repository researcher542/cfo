"""
Simulation for dual subgradient.
"""

"""
ds_flag: check the middle solution
"""
function ds_simulate(prob::AbsCfoProb, primal::DsPrimal; check_flag::Bool = true, debug::Bool = false, ds_flag::Bool = false, lb::Float64 = 0.0, predict_mode::AbstractPredictionMode, check_speed::Bool = false)
    (;net, ev, β0) = prob
    # cs_vec = net.cs_vec
    sol_vec = primal.sub_sol_vec
    obj::Float64 = 0.0
    beta_hat_vec::Vector{Float64} = []
    beta_vec::Vector{Float64} = [β0]
    tau_vec::Vector{Float64} = [0.0]
    beta_cs_vec::Vector{Float64} = [β0]
    tau_cs_vec::Vector{Float64} = [0.0]
    ci_vec::Vector{Float64} = []
    infeasible_flag::Bool = false
    obj1::Float64 = 0.0
    e_cost::Float64 = 0.0 # The charged energy 
    wasted_energy::Float64 = 0.0 ## The wasted_energy due to the battery capacity.
    for (isol,sol) in enumerate(sol_vec)
        (;path, t_vec, tc, tw) = sol
        path::Vector{Int} = sol.path
        t_vec::Vector{Float64} = sol.t_vec
        cost_v = (energy_cost_on_road(net, path[i], path[i+1], ev, t_vec[i]) for i in 1:length(path)-1)
        if check_speed
            for i in 1:length(path)-1
                min_t, max_t = get_minmax_t(net, path[i], path[i+1], ev)
                @show t_vec[i] min_t max_t
                if t_vec[i] < min_t || t_vec[i] > max_t
                    @warn "The time vector is not in the range of min and max time" isol i t_vec[i] min_t max_t
                end
            end
        end
        # @show isol, cost_v/3.6e6
        Δβ = charge_function(beta_vec[end], tc, ev.cap) - beta_vec[end]
        if isol > 1
            obj1 = cfo_objective(prob, path[1], beta_vec[end], tc, tau_vec[end]+tw, predict_mode)
            ci = obj1 / (Δβ + 1e-6)
            push!(ci_vec, ci)
            obj += obj1
            if debug 
                @show obj1 / 3.6e6
                @show tc / 60.0
                @show tw / 60.0
                @show Δβ / 3.6e6
                @show beta_vec[end] / 3.6e6
                @show tau_vec[end] / 60.0
                println("*"^80)
            end
        end
        # We set tolerance to 1min
        if (tau_vec[end] - sol.τ) / 60.0 > 1.0
            ds_flag && check_flag && @warn "scheduled larger actual time" isol tau_vec[end]/60.0 sol.τ/60.0
        end
        push!(tau_cs_vec, tau_vec[end])
        # The tolerance is 1 kWh
        if (beta_vec[end] - sol.β)/3.6e6 < -1.0
            ds_flag && check_flag && @warn "scheduled lower soc" isol beta_vec[end]/3.6e6 sol.β/3.6e6
        end
        push!(beta_cs_vec, beta_vec[end])

        Δτ = tw + tc
        e_cost += Δβ

        if debug
            @debug "At i=$isol, tau=$(tau_vec[end]), tw=$tw obj=$obj1 charged=$(Δβ) cs_idx=$(path[1])"
        end

        beta_hat_vec = vcat(beta_hat_vec, beta_vec[end] + Δβ)
        beta_tmp::Float64 = beta_vec[end] + Δβ
        for c1 in cost_v
            beta_tmp -= c1
            # @show c1
            beta_tmp = min(beta_tmp, ev.cap)
            # beta_vec = vcat(beta_vec, beta_tmp) 
            beta_vec = push!(beta_vec, beta_tmp) 
            # @show beta_vec[end] beta_tmp
        end
        ### Note: in this part, we need to simulate the beta step by step and consider the clamp from ev.cap
        # beta_vec = vcat(beta_vec, beta_vec[end] + Δβ .- cumsum(cost_v) )
        tau_vec = vcat(tau_vec, tau_vec[end] + Δτ .+ cumsum(t_vec) )
    end
    min_beta = minimum(beta_vec)
    if min_beta < -3.6e6 # less than 1 kWh
        check_flag && @warn (@sprintf "Infeasible with min_β=%.3f %%" min_beta/ev.cap*100) prob.src prob.des prob.T
        infeasible_flag = true
    end
    if tau_vec[end] > prob.T + 60.0 # more than 1 minute delay.
        check_flag && @warn (@sprintf "Infeasible with tau=%.2f > T=%.2f " tau_vec[end]/3600 prob.T/3600) prob.src prob.des prob.T
        infeasible_flag = true
    end

    ### account for the initial battery
    used_init_beta = β0 - beta_vec[end]
    e_cost += used_init_beta
    init_ci = get_init_battery_ci(prob, predict_mode)
    obj += init_ci * used_init_beta
    # @show used_init_beta / 3.6e6  beta_vec[end] / 3.6e6 init_ci
    
    # change it to kg, kWh, and minutes
    res = (obj=obj/3.6e6, beta_vec=beta_vec/3.6e6, tau_vec=tau_vec/60.0, infeasible_flag=infeasible_flag, e_cost=e_cost, tau_cs_vec=tau_cs_vec, beta_cs_vec=beta_cs_vec, lb=lb, beta_hat_vec=beta_hat_vec/3.6e6, ci_vec=ci_vec)
    return res
end


function get_init_battery_ci(prob::AbsCfoProb, predict_mode::AbstractPredictionMode)
    ta = get_price_ta(prob, prob.src, predict_mode)
    st = prob.start_time - Day(1)
    ci_vec = (get_price(ta, st, i*3600.0)[1] for i in 1:24)
    init_ci = mean(ci_vec)
    # @show ci_vec init_ci
    # init_ci = ci_vec[1]
    return init_ci
end