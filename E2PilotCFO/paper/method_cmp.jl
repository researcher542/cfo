"""
The code to compare the runtime of different methods.

We implement two methods: MIP formulation with Gurobi solver and the time-battery expanded network (TBExNet).

To run the code, you need to have the Gurobi solver installed and the correct lisence configured in your environment.
"""

import E2PilotCFO as cfo
import E2PilotCore as ep
using Dates, Logging, LoggingExtras, JLD2, Infiltrator
using Distributed, ArgParse

if Base.find_package("Gurobi") !== nothing
    using Gurobi
else
    @warn "Gurobi is not installed, so we cannot use the Gurobi solver."
end

cfo_solver = Base.get_extension(cfo, :SolverExt) # load the extension

cfo.init()


"""
"""
function reduced2full_idx(idx, prob, prob_us)
    n1 = ep.getnode(prob.net, idx)
    # new_idx = ep.closest_node_in_max_com(prob_us.net, n1)
    new_idx = ep.closest_node(prob_us.net, n1)
    idx2 = cfo.closest_cs_node(prob_us.net, new_idx)
    n2 = ep.getnode(prob_us.net, idx2)
    dis = ep.distance2d(n1, n2) 
    if dis > 0.1
        @warn "The distance between the two nodes is too large: $dis" idx idx2
        @exfiltrate
        @assert false
        return -1
    end
    return idx2
end


function get_one_region_cmp(region)
    # if region in ["486-485"]
    #     @warn "skip region $region. It might due to it is a hard case."
    #     return
    # end
    args = parse_arg()
    @info "================= $region ===================="
    if region == "test4node"
        prob = cfo.get_test_prob_small(region)
        prob_us = copy(prob)
    else
        prob = cfo.get_test_prob_us_reduced(region)
        prob_us = cfo.get_test_prob("oldmap")
    end

    prob_us.src = reduced2full_idx(prob.src, prob, prob_us)
    prob_us.des = reduced2full_idx(prob.des, prob, prob_us)
    # @info "For the od on prob_us: src: $(prob_us.src), des: $(prob_us.des)"
    # if prob_us.src == -1 || prob_us.des == -1
    #     @warn "The src or des is not in the reduced network."
    #     return
    # end
    if ep.nv(prob.net) > 10_000
        @warn "The network is too large, skip the region $region."
        return
    end

    # od_vec = cfo.read_od_data_vec(prob)
    verbose = 0
    (;ev) = prob_us
    g_alpha = 0.05
    maxtime = 60.0 * 60.0 ## change the maxtime to 20 minutes, for fast generation of results... We might need to change it later.
    mem_limit = 1024.0 * 16.0

    

    ds_option = cfo.DsOption(;N=-1, maxiter=100, maxtime=maxtime, β_lb=g_alpha*ev.cap, trace=true, rho=0e-1, debug=false, pos_cost=true, ds_net_type = cfo.DsNetwork(), early_break=false, verbose=verbose) # rho 


    alg_name0 = "fast-ds"
    fast_res = cfo.get_one_result_cmp(prob, prob_us, ds_option, alg_name0; ow_flag=false)

    if args["tbexnet"]
        alg_name0 = "tbexnet"
        T = fast_res.time * 1.2
        eps_t = min(T / 200.0, 60.0)
        tb_option = cfo.TbOption(;
            # eps_t = T / 200.0,
            eps_t = eps_t,
            eps_b = prob.ev.cap / 100.0,
            mem_limit = mem_limit,
        )
        tbexnet_res = cfo.get_one_result_cmp(prob, prob_us, tb_option, alg_name0; ow_flag=false, fast_res=fast_res)
    end

    if args["dsc"]
        alg_name0 = "ds-c"
        dsc_res = cfo.get_one_result_cmp(prob, prob_us, ds_option, alg_name0; ow_flag=false, fast_res=fast_res)
    end

    if args["grb"]
        if isnothing(cfo_solver)
            @warn "Gurobi is not installed, so we cannot use the Gurobi solver."
        else
            mip_option = cfo_solver.MIPOption(;
                verbose=1, 
                apx_type=[:ApxPiecewiseLinear, :ApxPolynomial][1],
                time_limit=maxtime,
                mem_limit=mem_limit,
            )
            alg_name0 = "grb"
            grb_res = cfo.get_one_result_cmp(prob, prob_us, mip_option, alg_name0; ow_flag=false, fast_res=fast_res)
        end
    end
   
    return
end

function parse_arg()

    s = ArgParseSettings()
    @add_arg_table s begin
        "--dsc"
        help = "If we should use the dsc method"
        arg_type = Bool
        default = true

        "--tbexnet"
        help = "If we should use the tbexnet method"
        arg_type = Bool
        default = true

        "--grb"
        help = "If we should use the grb method"
        arg_type = Bool
        default = true

        "--oneregion"
        help = "If we only want to run one region"
        arg_type = Bool
        default = false

        "--idx"
        help = "The index of the region to run"
        arg_type = Int
        default = 1

    end

    parsed_args = parse_args(ARGS, s)
    @show parsed_args

    return parsed_args
end

prob = cfo.get_test_prob("oldmap"); od_region_vec = cfo.get_cmp_od_region_vec(prob) 
region = od_region_vec[end]
# get_one_region_cmp(region)

log_prefix = "method-cmp-$(now())"
logger = ep.get_logger(log_prefix; console_level = Logging.Info)
args = parse_arg()

with_logger(logger) do
    @info "[$(now())] start getting results for method cmp."
    @info args
    # get_one_region_cmp("test4node")
    if args["oneregion"]
        region = od_region_vec[args["idx"]]
        stat = @timed get_one_region_cmp(region)
    else
        ## If we want to run single region for testing from Julia REPL.
        for region in ["171-179"]
            @info "start getting results for region $region."
            stat = @timed get_one_region_cmp(region)
            @info "get results for region $region, time: $(stat.time) seconds."
        end
    end
    
end
# exit(0)


# idx = cfo.find_od_idx(od_vec, prob.net.region)

# alg_name0 = "tbexnet"
# cfo_res = cfo.get_one_result_cmp(prob, option, alg_name0; ow_flag=true)






