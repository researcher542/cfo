import E2PilotCore as ep
import E2PilotCFO as cfo
import E2PilotCFO
import E2PilotCore
using LoggingExtras
using Dates
using TimerOutputs
using Geodesy
using Printf
using CSV, DataFrames, JLD2

cfo.init()

# ENV["JULIA_DEBUG"] = "Main,E2PilotCore,E2PilotCFO"
ENV["JULIA_DEBUG"] = "Main,E2PilotCore"
sf = ep.safehouse
reset_timer!(ep.g_to)

# If true, 
cs_nei_flag = false

sep_cs_node_flag = false
# if true, we will separate the charging station node to a different node and add it to the graph.
if !isdefined(Main, :sep_cs_node_flag)
    sep_cs_node_flag = true
end
