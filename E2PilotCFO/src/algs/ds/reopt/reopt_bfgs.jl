
import Optim 

"""
Re-optimize primal solution with BFGS.
!!!Do not work well
"""
function ds_reopt_bfgs(prob::AbsCfoProb, primal::DsPrimal, option::DsOption)

    lb, ub = ds_reopt_cs_net_get_lbub(prob, primal, option)    
    x0 = ds_reopt_cs_net_primal2tvec(primal)
    x0 = clamp.(x0, lb, ub)
    # x1 = (lb+ub)./2

    function f(x)
        ret = ds_reopt_cs_net_objective(x, prob, primal, option)
        # @show x ret
        return ret
    end
    function g!(G, x)
        ret = FiniteDiff.finite_difference_gradient!(G, f, x)
        # @warn "" x ret
        # for i in eachindex(G)
        #     G[i] = ret[i]
        # end
        return ret
    end
    inner_opt = Optim.BFGS()
    fmin_opt = Optim.Fminbox(inner_opt)
    opt_option = Optim.Options(
        show_trace=false
    )
    res = Optim.optimize(f, g!, lb, ub, x0, fmin_opt, opt_option)
    # res = Optim.optimize(f, g!, x0, inner_opt, opt_option)
    cost = Optim.minimum(res)
    x = Optim.minimizer(res)

    primal_new = ds_reopt_cs_net_tvec2primal(x, primal, prob, option)
    return cost, primal_new
end