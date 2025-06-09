


"""
Re-optimize primal with black-box optiomization. 
The travelling time and cost relies on the the cs_net. 
Hope that we can speed up the optimization.
"""
function ds_reopt_cs_net_bbo(prob::AbsCfoProb, primal::DsPrimal, option::DsOption)
    # (;net, ev) = prob

    lb, ub = ds_reopt_cs_net_get_lbub(prob, primal, option)    

    x0 = ds_reopt_cs_net_primal2tvec(primal)
    x0 = clamp.(x0, lb, ub)
    x1 = (lb+ub)./2
    # @show lb ub x0

    method = SA[:generating_set_search, :adaptive_de_rand_1_bin_radiuslimited, :adaptive_de_rand_1_bin][2]

    function f(x)
        obj = ds_reopt_cs_net_objective(x, prob, primal, option)
        return obj
    end
    @timeit g_to "reopt_cs_net" res = BBO.bboptimize(f, [x0,x1]; SearchRange=[(l,u) for (l,u) in zip(lb, ub)], Method = method, TraceMode= :silent, MaxFuncEvals= 5_000_00) #, MaxTime = 10)

    cost = BBO.best_fitness(res)
    x = BBO.best_candidate(res)

    primal_new = ds_reopt_cs_net_tvec2primal(x, primal, prob, option)

    return cost, primal_new
end
