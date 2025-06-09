
import QuadDIRECT as qd

"""
Solve the non-convex sub-problem to get σ_v^i with QuadDIRECT
i: i-th charging stop during the trip
return (cost, (tc,tw,β,τ))
"""
function ds_get_cs_cost_direct(prob::AbsCfoProb, primal::DsPrimal, dual::DsDual, istage::Int, cs_road_idx::Int, ta::TS.TimeArray, option::DsOption)
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

    function f(x)  
        # obj = ds_sigma_objective(x, primal, dual, i, ta, prob.start_time, B, option)
        if x[1] == 0.0
            x[2] = 0.0
        end
        obj = ds_sigma_objective(prob, x, primal, dual, istage, cs_road_idx, option)
        if isnan(obj) || isinf(obj)
            return 1e20
        end
        return obj
    end

    splits = [[lb[i], (lb[i]+ub[i])/2, ub[i]] for i in 1:4] 
    # splits = [zeros(nsplits) for i in 1:4]
    # for isplit in 1:nsplits
    #     for idx in 1:4
    #         flag = (isplit & (2^(idx-1))) != 0
    #         if flag
    #             splits[idx][isplit] = lb[idx]
    #         else
    #             splits[idx][isplit] = ub[idx]
    #         end
    #     end
    # end

    # splits = (lb, ub, x0)
    # splits = x0
    minwidth = [1.0, 1.0, B/1e6, 5.0]

    print_interval = [1_000_0, typemax(Int)][2]
    # maxevals = 10_000
    # maxevals = 2500
    maxevals = 1000
    max_tol_counter = 10
    # nquasinewton = 1000
    root, x0 = qd.analyze(f, splits, lb, ub; 
        minwidth=minwidth, 
        print_interval=print_interval, maxevals = maxevals, 
        rtol=1e-3, max_tol_counter=max_tol_counter, nquasinewton=maxevals*2)
    box = qd.minimum(root)

    x = qd.position(box, x0)
    if x[1] == 0.0
        x[2] = 0.0
    end
    cost = f(x)

    # @show lb1 ub1 prob.objtype x cost f([0.0, x[2], x[3], x[4]]) f([x[1], x[2], x[3], x[4]])
   

    return (cost, x)
end