"""
Get the all the results (var. ddl., var. year...) for one OD
"""

include("testdata/pre.jl")
using Dates, ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--region"
    help = "the region to test"
    default = "oldmap"

    "--idx"
    help = "the index of OD pair"
    arg_type = Int
    default = 1

    "--ow"
    help = "overwrite the existing result"
    arg_type = Bool
    default = false

    "--dry_run"
    help = "dry run"
    arg_type = Bool
    default = false

    "--run_alpha"
    help = "run alpha related experiments"
    arg_type = Bool
    default = true

    "--run_t_ratio"
    help = "run t_ratio related experiments"
    arg_type = Bool
    default = true

    "--run_future"
    help = "run future related experiments"
    arg_type = Bool
    default = true

    "--run_noise"
    help = "run ci noise related experiments"
    arg_type = Bool
    default = true

    "--base_only"
    help = "Run the basic cases only, this will disable all other flags"
    arg_type = Bool
    default = false
end

    # Parse the arguments
parsed_args = parse_args(ARGS, s)

region = parsed_args["region"]


ow_flag = parsed_args["ow"]
dry_run = parsed_args["dry_run"]
run_alpha = parsed_args["run_alpha"]
run_t_ratio = parsed_args["run_t_ratio"]
run_future = parsed_args["run_future"]
run_noise = parsed_args["run_noise"]
idx::Int = parsed_args["idx"]
base_only = parsed_args["base_only"]

if base_only
    run_alpha = run_t_ratio = run_future = run_noise = false
end

@info "parsed_args:" idx ow_flag dry_run run_alpha run_t_ratio run_future run_noise base_only
verbose = 0


if length(ARGS) >= 1
    disable_timer!(cfo.g_to)
    ENV["JULIA_DEBUG"] = ""
    @show ARGS
else
    # disable_timer!(cfo.g_to)
    # if not called from cmd, then it must be testing enviroment.
    ENV["JULIA_DEBUG"] = "Main,E2PilotCore,E2PilotCFO"
    dry_run = true
    ow_flag = false
    idx = 1
    verbose = 1
    run_alpha = false
    run_t_ratio = false
    run_future = false
    run_noise = true
    region = ["oldmap", "eu"][1]
end

if !isdefined(Main, :prob)  || isnothing(prob)
    prob = cfo.get_test_prob(region)
end
(; ev, net) = prob
g_alpha = 0.05
option = cfo.DsOption(;N=-1, maxiter=100, maxtime=3600.0, β_lb=g_alpha*ev.cap, trace=true, rho=0e-1, debug=false, pos_cost=true, ds_net_type = cfo.DsNetwork(), early_break=false, verbose=verbose) # rho candidate: 0.5 or 1.0
β0 = ev.cap
ice_veh = ep.get_veh()

@show idx, ow_flag, gethostname()
flush(stdout)
# ow_flag = true

log_prefix = "idx=$idx.oneOD.$(region).-$(now())"
logger = ep.get_logger(log_prefix; console_level = Logging.Info)

# prob = cfo.CfoProb(;fix_data=Ref(fix_data), src=src, des=des, β0=β0, T=3600.0 * 100.0, start_time=start_time)
with_logger(logger) do
    @info "[$(now())] start getting one src des idx=$idx."
    @info "parsed_args:" idx ow_flag dry_run run_alpha run_t_ratio run_future run_noise base_only
    global res_vec = cfo.get_one_src_des(prob, idx, ice_veh, option; ow_flag=ow_flag, dry_run=dry_run, 
        run_alpha=run_alpha, 
        run_t_ratio=run_t_ratio, 
        run_future=run_future,
        run_noise=run_noise
        )
    if dry_run
        for arg in res_vec
            alg_name = arg[2]
            @show alg_name
        end
    end
    @info "[$(now())] idx=$idx done."
end
@info "[$(now())] done."