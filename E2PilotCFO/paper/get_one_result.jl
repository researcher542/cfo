"""
test for get one result
"""

# include("testdata/largedata.jl") ## this is the map used in the paper.
# include("testdata/TXOH.jl") ## use smaller data set for testing.
# ENV["JULIA_DEBUG"] = "Main,E2PilotCore,E2PilotCFO"

include("testdata/pre.jl")

import E2PilotCFO as cfo
import E2PilotCore as ep
using Dates

cfo.init()

# prob = cfo.get_test_prob("oldmap") ## the complete US network.
prob = cfo.get_test_prob("TXOH")
(; ev, net, start_time) = prob

g_alpha = 0.05
option = cfo.DsOption(;N=-1, maxiter=100, maxtime=3600.0, β_lb=g_alpha*ev.cap, trace=true, rho=0e-1, debug=false, pos_cost=true, ds_net_type = cfo.DsNetwork(), early_break=false, verbose=0) 
ice_veh = ep.get_veh()


if net.region == "oldmap" ## if it is the large map.
    alpha = 0.05
    t_ratio = 1.20
    idx = 1
    src_des_vec = ep.read_faf_data_net(net, 1000)
    src, des = src_des_vec[idx]
    src, des = cfo.closest_cs_node.((net,), (src,des))
else
    alpha = 0.05
    t_ratio = 1.20
    idx = 1
    src_des_vec = ep.read_faf_data_net(net, 1000)
    src, des = src_des_vec[idx]
    src, des = cfo.closest_cs_node.((net,), (src,des))
end
# src, des = 17745, 20970
## get the fast-ds result.
alg_name = "fast-ds"
fast_res = nothing
res = cfo.CfoResult(; alg="", src=src, des=des, idx=idx, start_time=start_time, cap=ev.cap, β0=prob.β0)
res.region = net.region

args = (prob, alg_name, deepcopy(res), t_ratio, option, fast_res, ice_veh, idx, false)
fast_ds_res = cfo._get_one_result(args)

show(ep.g_to)
println()

## get for ds-c

now_str = Dates.format(now(), "yyyy-mm-dd-HH.MM.SS")
logger = ep.get_logger("get_one_result-$(now_str)")

with_logger(logger) do
    alg = "ds-c"
    prob.objtype = [cfo.ObjCarbon(), cfo.ObjEnergy()][1]
    alg_name = cfo.get_result_key(alg, t_ratio, alpha, prob.start_time)
    option.β_lb = alpha * ev.cap
    res = cfo.CfoResult(; alg="", src=src, des=des, idx=idx, start_time=start_time, cap=ev.cap, β0=prob.β0)
    res.region = net.region
    args = (prob, alg_name, deepcopy(res), t_ratio, option, fast_ds_res, ice_veh, idx, true)
    @time global ds_res = cfo._get_one_result(args)
    # @show ds_res.cost
end

show(ep.g_to)
println()

# results for the
println("*"^40)
println("Results for the carbon-optmized solution.")
println("Carbon footprint: $(ds_res.cost / 1e6) kg")
println("Time: $(ds_res.time / 3600) hours")
println("Energy: $(ds_res.e_cost / 3.6e6) kWh")

## results for the fastest path
println("*"^40)
println("Results for the fastest path.")
println("Carbon footprint: $(fast_ds_res.cost / 1e6) kg")
println("Time: $(fast_ds_res.time / 3600) hours")
println("Energy: $(fast_ds_res.e_cost / 3.6e6) kWh")


# args = (prob, alg_name, deepcopy(res), g_t_ratio, deepcopy(option), fast_ds_res, ice_veh, idx)
