
import QuadDIRECT as qd 

"""
Re-optimize primal solution with QuadDIRECT Library
"""
function ds_reopt_direct(prob::AbsCfoProb, primal::DsPrimal, option::DsOption)

    lb, ub = ds_reopt_cs_net_get_lbub(prob, primal, option)    
    nx = length(lb)
    ub = ub .+ 1e-6 * ones(nx)
    # x0 = ds_reopt_cs_net_primal2tvec(primal)
    # x0 = clamp.(x0, lb, ub)
    x0 = (lb+ub)./2

    function f(x)
        # ret = ds_reopt_cs_net_objective(x, prob, primal, option)
        obj, b_cons_vec, t_cons = ds_reopt_cs_net_obj_cons(x, prob, primal, option)
        rho = 1e8
        b_cons0 = rho * sum((pos(b_cons_vec)).^2)
        t_cons0 = rho * pos(t_cons)^2
        ret = obj + b_cons0 + t_cons0 
        return ret
    end

    splits = [[lb[i], (lb[i]+ub[i])/2, ub[i]] for i in 1:nx] 
    minwidth = ones(nx)
    print_interval = [1_000_0, typemax(Int)][2]
    # maxevals = 10_000
    maxevals = 20_000
    max_tol_counter = 50
    # nquasinewton = 1000
    root, x0 = qd.analyze(f, splits, lb, ub; minwidth=minwidth, print_interval=print_interval, maxevals = maxevals, rtol=1e-3, max_tol_counter=max_tol_counter)
    box = qd.minimum(root)
    x = qd.position(box, x0)
    cost = f(x)
   
    # @show ret = ds_reopt_cs_net_obj_cons(x, prob, primal, option, true)
    # @timeit g_to "tvec2primal" 
    primal_new = ds_reopt_cs_net_tvec2primal(x, primal, prob, option)
    
    return cost, primal_new
end