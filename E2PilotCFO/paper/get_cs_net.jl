

using LoggingExtras
include("testdata/pre.jl")

if Sys.isapple() || occursin("sujy", homedir()) || occursin("ubuntu", homedir())
    @info "getting cs_net in local computers."
    # prob = cfo.get_test_prob_small("test4node")
    # prob = cfo.get_test_prob_small("testsmall")
    # prob = cfo.get_test_prob("TXOH")
    prob = cfo.get_test_prob("eu")
    # include("testdata/mideast.jl")
    # include("testdata/test4.jl")
    # include("testdata/test6.jl")
    # include("testdata/TXOH.jl")
    # include("testdata/apx.jl")
    # include("testdata/AZWA.jl")
    ow_flag = false
else
    @info "getting cs_net in remote server."
    # prob = cfo.get_test_prob_large()
    prob = cfo.get_test_prob_small("test4node"; B = 3.6e9)
    # include("testdata/largedata.jl")
    # include("testdata/mideast.jl")
    # include("testdata/TXOH.jl")
    # include("testdata/AZWA.jl")
    # include("testdata/test4.jl")
    ow_flag = false
end

cs_nei_flag = false
prob.ev.cap = 3.6e9
# cfo.set_no_regen!(ev) ## 

# T = 12*3600.0
# path = ep.shortest_path(ep.Astar(), cs_net, src, des)
# step0 = ep.paso(net, src, des, Inf; veh=ev, visitonce=true)
logger = ep.get_logger("get_cs_net-$(now())")
# @assert false

with_logger(logger) do 
    global cs_net
    @info "getting cs_net for region $(prob.net.region) with nthread=$(Threads.nthreads()) with cs_nei_flag=$cs_nei_flag"
    @timeit ep.g_to "get_cs_net" cs_net = cfo.compute_cs_net(prob.net, prob.ev; ow_flag=ow_flag, cs_nei_flag=cs_nei_flag)
    cfo.save_cs_net(prob.net.region, cs_net, cs_nei_flag; tmp=false)
end

show(ep.g_to)
@info "done"