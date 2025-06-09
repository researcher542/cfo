
# struct to determine use which network type
abstract type AbsDsNetType end
struct DsNetwork <: AbsDsNetType end
struct DsCsNet <: AbsDsNetType end

@with_kw mutable struct DsOption{NT <: AbsDsNetType}
    ds_net_type::NT = DsNetwork()  # struct to determine use which network type to use 
    N::Int = -1 # number of charging stations 
    maxiter::Int = 100
    maxtime::Float64 = 1000.0 # max time in seconds
    β_lb::Float64 = 0.0 # The β lower bound when entering the charging station. used to eliminate the effect of regenerative system.
    trace::Bool = false
    pos_cost::Bool = false # enforce positive cost to avoid negtive cycle. This may lead to sub-optimality.

    primal0::DsPrimal = DsPrimal() # The initial primal for warm-start

    debug::Bool = false
    early_break::Bool = true

    # The multipilers to ensure numerical stability and scale the problem well
    # Ideally, we consider 1 min <-> 1 kWh
    b_mul::Float64 = 1/3.6e6 # from J change to kWh
    t_mul::Float64 = 1/60.0 # from second change to minutes
    rho::Float64 = 0.0 # The panalty of "smooth" the iteration.

    # The preprocessed information
    fast_t_vec::Vector{Float64} = [] # the vector that store the fastest time from src to each node
    fast_t_rev_vec::Vector{Float64} = [] # the vector that store the fastest time from each node to des, rev means reverse

    verbose::Int = 1
    kpath_flag::Bool = true # If we want to use kpath to generate a feasible path at the begining.
end