
"""
Solve the non-convex sub-problem to get σ_v^i with NLopt
i: i-th charging stop during the trip
return (cost, (tc,tw,β,τ))
"""
function ds_get_cs_cost_nlopt(prob::AbsCfoProb, primal::DsPrimal, dual::DsDual, istage::Int, cs_road_idx::Int, ta::TS.TimeArray, option::DsOption)
    B = prob.ev.cap
    # (;λ,μ) = dual

    lb1,ub1 = ds_get_cs_var_range(prob, primal, istage, cs_road_idx, option)
    lb = Vector(lb1)
    ub = Vector(ub1)
    x0 = (lb+ub)./2
    if lb[4] > ub[4]
        # @show istage,cs_road_idx,lb[4],ub[4]
        return (Inf, Tuple(lb))
    end
    @assert (all(lb .<= ub)) 

    NLopt.srand(0)

    function f(x, grad)  
        # obj = ds_sigma_objective(x, primal, dual, i, ta, prob.start_time, B, option)
        obj = ds_sigma_objective(prob, x, primal, dual, istage, cs_road_idx, option)
        return obj
    end

    function _f(x, grad)
        try f(x, grad) 
        catch e
            bt = catch_backtrace()
    	    showerror(stdout, e, bt)
    	    rethrow(e)
        end
    end


    flag = 1
    if flag == 1
        method = SA[:GN_DIRECT_L, :GN_AGS, :GN_ISRES, :GN_ESCH, :GN_CRS2_LM, :LN_COBYLA][1]
        opt = NLopt.Opt(method, length(x0))
        opt.lower_bounds = lb
        opt.upper_bounds = ub
        # opt.xtol_rel = 1e-8
        opt.min_objective = _f
        # opt.maxtime = 10.0
        opt.maxeval = 200_000
    elseif flag == 2

    end

    (cost,x,ret) = NLopt.optimize(opt, x0)
    # @show ret cost x
    if ret == :FORCED_STOP
        error("force stoped in nlopt")
    end

    return (cost, x)
end