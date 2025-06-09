
import NLopt 

"""
Re-optimize primal solution with NLOPT Library
NOTE: should only use direvative-free method as the direvative of our objective is not well-defined.
"""
function ds_reopt_nlopt(prob::AbsCfoProb, primal::DsPrimal, option::DsOption, method = nothing)

    lb, ub = ds_reopt_cs_net_get_lbub(prob, primal, option)    
    x0 = ds_reopt_cs_net_primal2tvec(primal)
    x0 = clamp.(x0, lb, ub)
    # x1 = (lb+ub)./2

    function f(x)
        # ret = ds_reopt_cs_net_objective(x, prob, primal, option)
        ret = ds_reopt_cs_net_obj_cons(x, prob, primal, option)[1]
        return ret
    end

    function myfunc(x, grad)
        if length(grad) > 0
            g!(grad, x)
        end
        obj = f(x) 
        return obj
    end

    function _myfunc(x, grad)
        try myfunc(x, grad)
        catch e
    	    bt = catch_backtrace()
    	    showerror(stdout, e, bt)
    	    rethrow(e)
        end
    end

    function mycons(result::Vector, x::Vector, grad::Matrix)
        obj, b_cons_vec, t_cons = ds_reopt_cs_net_obj_cons(x, prob, primal, option)
        n_b_cons = length(b_cons_vec)
        result[1:n_b_cons] = b_cons_vec
        result[n_b_cons+1] = t_cons
    end

    function _mycons(result, x, grad)
        try mycons(result, x, grad)
        catch e
    	    bt = catch_backtrace()
    	    showerror(stdout, e, bt)
    	    rethrow(e)
        end
    end

    if isnothing(method)
        # method = [:GN_DIRECT_L, :GN_AGS, :GN_ISRES, :GN_CRS2_LM, :GN_ESCH][4]
        # global with constraint
        method = [:GN_ISRES][1]
        # local optimizer
        # method = [:LN_COBYLA, :LN_BOBYQA][2]
    end


    NLopt.srand(0)
    maxtime = 100.0
    maxeval = 5_000_00

    if false
        opt = NLopt.Opt(:AUGLAG, length(x0))
        local_method = [:LN_COBYLA, :LN_BOBYQA, :LN_NELDERMEAD, :LN_SBPLX, :LN_NEWUOA, :LN_PRAXIS][4]
        # local_method = [:LN_SBPLX]
        local_opt = NLopt.Opt(local_method, length(x0))
        opt.local_optimizer = local_opt
        # local_opt.lower_bounds = lb
        # local_opt.upper_bounds = ub
        min_out_iter = 10
        # local_opt.maxeval = round(Int, maxeval / min_out_iter)
        # local_opt.maxtime = maxtime / min_out_iter
        # local_opt.xtol_rel = 1e-4
    else
        opt = NLopt.Opt(method, length(x0))
    end

    opt.lower_bounds = lb
    opt.upper_bounds = ub

    # opt.xtol_rel = 1e-8
    opt.min_objective = _myfunc
    N = get_N(primal)
    cons_tol_vec = zeros(N+3)
    NLopt.inequality_constraint!(opt, _mycons, cons_tol_vec)
    opt.maxeval = maxeval
    opt.maxtime = maxtime

    (cost,x,ret) = NLopt.optimize(opt, x0)
    # @warn "" prob.objtype cost 
    (option.verbose > 0) && @show ret
    if ret == :FORCED_STOP
        error("forced stopped.")
    end

    # ds_reopt_cs_net_objective(x, prob, primal, option, true)
   
    # @show ret = ds_reopt_cs_net_obj_cons(x, prob, primal, option, true)
    @timeit g_to "tvec2primal" primal_new = ds_reopt_cs_net_tvec2primal(x, primal, prob, option)
    
    return cost, primal_new
end