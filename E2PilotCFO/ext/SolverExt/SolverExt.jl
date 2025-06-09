
"""
Use existing MIP solvers to solve the CFO problem.
"""

module SolverExt

using JuMP
using Parameters, Graphs

import E2PilotCFO as cfo
import E2PilotCore as ep

import E2PilotCFO: g_min_charge_time, g_max_charge_time, g_min_wait_time, g_max_wait_time, charge_function, cf_inv, get_minmax_t, get_price, get_price_ta, DsSubsolution, DsPrimal

include("jump.jl")



function test_cfo()
    @info "This is a test message in solver CFO"
end

end