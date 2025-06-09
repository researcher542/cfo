
include("bbo_sigma.jl")
include("bnb_sigma.jl")
include("nlopt_sigma.jl")
include("direct_sigma.jl")

"""
Get the augmented cost for each node of charging station

return:    cs_cost_vec[istage][ics] = cost
"""
function ds_get_cs_cost_vec(prob::AbsCfoProb, cs_vec::Vector{ChargeStation}, primal::DsPrimal, dual::DsDual, option::DsOption)
    (;src, des, net) = prob
    n_cs = length(cs_vec)
    N = option.N
    cs_cost_vec = [zeros(n_cs) for _ in 1:N]
    cs_idx_vec = [ics for ics in eachindex(cs_vec)]
    function func(ics) 
        cs = cs_vec[ics]
        cs_road_idx = cs.idx
        price_ta = get_price_ta(prob, cs_road_idx, prob.predict_mode)
        cost_v::Vector{Float64} = [ds_get_cs_cost(prob, primal, dual, istage, cs_road_idx, price_ta, option)[1] for istage in 1:N]
        return cost_v
    end
    # task_vec = [Threads.@spawn func(ics) for ics in cs_idx_vec]
    # res_vec = fetch.(task_vec)
    res_vec = e2map(func, cs_idx_vec, false)
    for ics in cs_idx_vec
        for istage in 1:N
            cs_cost_vec[istage][ics] = res_vec[ics][istage]
        end
    end
    return cs_cost_vec
end


"""
We decide the range of decision variables based on previous decisions
"""
function ds_get_cs_var_range(prob::AbsCfoProb, primal::DsPrimal, istage::Int, cs_road_idx::Int, option::DsOption)

    (;beta_cs_vec, tau_cs_vec, sub_sol_vec) = primal
    (;T, ev) = prob
    (;β_lb, fast_t_vec, fast_t_rev_vec) = option
    # the index of charging station in the network
    # cs_road_idx = sub_sol_vec[i].path[end]
    B = ev.cap
    # ratio = option.cs_var_ratio

    beta_cs = beta_cs_vec[istage+1]
    tau_cs = tau_cs_vec[istage+1]
    last_tc = sub_sol_vec[istage].tc
    last_tw = sub_sol_vec[istage].tw

    # tc_lb = max(last_tc - g_max_charge_time*ratio, g_min_charge_time)
    # tc_ub = min(last_tc + g_max_charge_time*ratio, g_max_charge_time)
    tc_lb = g_min_charge_time
    tc_ub = g_max_charge_time
    tc_ub = min(tc_ub, cf_inv(β_lb, B, B))

    # beta_lb = max(beta_cs - B*ratio , option.β_lb)
    # beta_ub = min(beta_cs + B*ratio , B)
    beta_lb, beta_ub = option.β_lb, B

    if !isempty(fast_t_vec)
        min_tau = fast_t_vec[cs_road_idx]
        max_tau = T - fast_t_rev_vec[cs_road_idx]
    else
        min_tau, max_tau = 0.0, T
    end
    tau_lb, tau_ub = min_tau, max_tau

    tw_lb, tw_ub = g_min_wait_time, g_max_wait_time

    if prob.des == cs_road_idx
        tw_lb, tw_ub = 0.0, 1e-4
        tc_lb, tc_ub = 0.0, 1e-4
    end
    # tau_lb, tau_ub = 0.0, T
    # if min_tau <= max_tau
    #     tau_lb = clamp(tau_cs - T*ratio, min_tau, max_tau)
    #     tau_ub = clamp(tau_cs + T*ratio, min_tau, max_tau)
    # else
    #     tau_lb, tau_ub = min_tau, max_tau
    # end

    # @assert tau_lb <= tau_ub "ratio=$ratio  tau_cs=$tau_cs T=$T tau_lb=$tau_lb tau_ub=$tau_ub"

    # tc, tw, b, tau
    lb = MVector{4}([tc_lb, tw_lb, beta_lb, tau_lb])
    ub = MVector{4}([tc_ub, tw_ub, beta_ub, tau_ub])
    
    return lb,ub
end


"""
x: [tc, tw, β, τ]
p: the coefficnet for the above 
"""
function ds_sigma_objective(prob::AbsCfoProb, x::AbstractVector, primal::DsPrimal, dual::DsDual, istage::Int, road_idx::Int, option::DsOption)
    B = prob.ev.cap
    (tc::Float64, tw::Float64, β::Float64, τ::Float64) = x
    (;t_mul, b_mul, rho, N) = option
    (;λ,μ) = dual
    (;beta_cs_vec, tau_cs_vec) = primal
    # c_tc::Float64, c_tw::Float64, c_β::Float64, c_τ::Float64, c_Δb::Float64, start_time, B, ta = p
    Δβ = charge_function(β, tc, B) - β
    obj::Float64 = cfo_objective(prob, road_idx, β, tc, τ+tw, prob.predict_mode)
    # λ[i+1]
    λi1 = (istage == N+1) ? 0.0 : λ[istage+1]
    μi1 = (istage == N+1) ? 0.0 : μ[istage+1]
    c_obj = obj * b_mul
    if prob.objtype == ObjTime()
        c_obj = (tc+tw) * t_mul
    end
    c_tcw = λi1 * (tw+tc) * t_mul
    c_tau = (λi1 - λ[istage]) * τ * t_mul
    c_beta = (μ[istage]-μi1)*β* b_mul
    c_del_beta = -μi1 * Δβ * b_mul

    last_beta = primal.sub_sol_vec[istage].β
    last_tau = primal.sub_sol_vec[istage].τ

    # beta_cs = max(beta_cs_vec[i+1], 0.0)
    beta_cs = beta_cs_vec[istage+1]
    tau_cs = tau_cs_vec[istage+1]
    last_tc = primal.sub_sol_vec[istage].tc
    last_tw = primal.sub_sol_vec[istage].tw
    c6 = rho * b_mul^2 * ( beta_cs - β )^2
    c7 = rho * t_mul^2 * ( tau_cs - τ )^2
    c8 = rho * t_mul^2 * ( last_tc - tc )^2
    c9 = rho * t_mul^2 * ( last_tw - tw )^2

    obj = c_obj + c_tcw + c_tau + c_beta + c_del_beta + c6 + c7 + c8 + c9
    if option.debug 
        @show c_obj c_tcw c_tau c_beta c_del_beta obj Δβ; println()
    end
    return obj
end

# ds_get_cs_cost(prob, dual, i, ta, option) = ds_get_cs_cost_bbo(prob, dual, i, ta, option)
function ds_get_cs_cost(prob::AbsCfoProb, primal::DsPrimal, dual::DsDual, istage::Int, cs_road_idx::Int, ta::TS.TimeArray, option::DsOption)      
    # @time "nlopt" 
    # (cost_nlopt, x_nlopt) = ds_get_cs_cost_nlopt(prob, primal, dual, istage, cs_road_idx, ta, option)
    # @time "direct" 
    (cost_d, x_d) = ds_get_cs_cost_direct(prob, primal, dual, istage, cs_road_idx, ta, option)
    # rel_err = (cost_d - cost_nlopt)/cost_nlopt
    # @debug rel_err cost_d cost_nlopt 
    # @show cost_nlopt
    # @time "bbo" (cost_bbo, x_bbo) = ds_get_cs_cost_bbo(prob, primal, dual, istage, cs_road_idx, ta, option)
    # rel_err = (cost_nlopt-cost_bbo)/max(cost_nlopt, cost_bbo)
    # ((abs(rel_err ) > 1e-3) && abs(cost_d) > 1e-6) &&  @show "-----" cost_nlopt cost_d rel_err "-----"
   
    (cost, x) = (cost_d, x_d)
    # (cost, x) = (cost_nlopt, x_nlopt)
    # return cost_bbo, x_bbo
    # return cost_nlopt, x_nlopt
    return cost, x
end



