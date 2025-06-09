
"""
Solve the non-convex sub-problem to get σ_v^i with BBO
i: i-th charging stop during the trip
return (cost, (tc,tw,β,τ))
"""
function ds_get_cs_cost_bbo(prob::AbsCfoProb, primal::DsPrimal, dual::DsDual, i::Int, cs_road_idx::Int, ta::TS.TimeArray, option::DsOption)
    B = prob.ev.cap
    (;λ,μ) = dual
    # T = prob.T
    # N = length(λ) - 1

    lb1,ub1 = ds_get_cs_var_range(prob, primal, i, cs_road_idx, option)
    lb = Vector(lb1)
    ub = Vector(ub1)
    
    x0 = (lb+ub)./2

    # directly use BBO, maxTime might be too small...
    method = SA[:generating_set_search, :adaptive_de_rand_1_bin_radiuslimited][1]
    # f = x->ds_sigma_objective(x, primal, dual, i, ta, prob.start_time, B, option)
    f = x->ds_sigma_objective(prob, x, primal, dual, i, cs_road_idx, option)
    res = BBO.bboptimize(f, x0; SearchRange=[(l,u) for (l,u) in zip(lb, ub)], Method = method, TraceMode= :silent , MaxTime = 1000e-3)

    cost = BBO.best_fitness(res)
    x = BBO.best_candidate(res)
    # error("")

    # @show cost x p[1:5]
    return (cost, x)
end