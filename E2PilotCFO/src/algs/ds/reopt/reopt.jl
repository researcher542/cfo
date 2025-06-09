
include("reopt_utils.jl")
include("reopt_obj_cons.jl")
include("reopt_bbo.jl")
include("reopt_bbo_cs_net.jl")
include("reopt_bfgs.jl")
include("reopt_nlopt.jl")
include("reopt_dp.jl")
include("reopt_direct.jl")

include("update_primal.jl")

function print_timed_stat(stat)
    gcstat = stat.gcstats
    for sym in propertynames(gcstat)
        @show sym, getproperty(gcstat,sym)
    end
    @show stat.time
    @show stat.bytes/1e6
    @show stat.gctime
    @show stat.compile_time
    @show stat.recompile_time
end

"""
re-optimize the primal to get a feasible solution. This is the high-level function to be called.

return: obj, primal
"""
function ds_reopt(prob::AbsCfoProb, primal0::DsPrimal, option::DsOption) 
    # if !ds_check_cs_path_feasible(prob, primal)
        # @warn "The cs path in the primal is not feasible!" get_cs_path(primal)
        # return Inf, primal
    # end
    primal = ds_make_cs_path_feasible(prob, primal0)
    @assert(ds_check_cs_path_feasible(prob, primal))
    if option.N != -1 && get_N(primal) > option.N 
        if option.verbose > 0
            # @warn "The number of charging stops is larger than option.N=$(option.N), return Inf." get_N(primal)
        end
        return Inf, primal
    end

    (obj_tmp, primal_tmp) = try
        stat = @timed obj_tmp, primal_tmp = ds_reopt_dp(prob, primal, option)
        if option.verbose > 0
            # print_timed_stat(stat) 
        end
        (obj_tmp, primal_tmp)
    catch e
        @warn "Error happens during ds_reopt, return " e prob primal0 option prob.src prob.des prob.T
        rethrow(e)
        return Inf, primal
    end
    
    debug = false
    if isinf(obj_tmp) && debug
        @timeit "ds_reopt_nlopt" obj_nlopt, primal_nlopt = ds_reopt_nlopt(prob, primal, option)
        @warn "reopt_dp failed, try to reopt with nlopt, got nlopt_obj=$(obj_nlopt), dp_obj=$(obj_tmp)"
        return obj_nlopt, primal_nlopt
    else
        return obj_tmp, primal_tmp
    end
end