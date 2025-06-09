
@with_kw mutable struct PbmState <: AbsDualState
    k::Int = 1
   
    # debug information
    cons_vio::Float64 = 0.0 # The constraint violation
    obj::Float64 = 0.0 # The raw objective
    lag_obj::Float64 = 0.0 # the Lagrangian objective
    dual::DsDual = DsDual(0)
    primal::DsPrimal = DsPrimal()

    # _c stands for the current best Lagrangian objective
    theta0::Float64 = - Inf
    theta_c::Float64 = -Inf # The current best Lagrangian objective
    dual_c::DsDual = DsDual(0) # current best dual
    primal_c::DsPrimal = DsPrimal() # current best primal

    # _f stands for the current best feasible point
    cons_vio_f::Float64 = Inf # The current best constraint violation.
    obj_f::Float64 = Inf # The current best feasible objective
    obj_f0::Float64 = Inf # The initial feasible objective by reoptimization
    theta_f::Float64 = Inf # The current best feasible Lagrangian
    # dual_f::DsDual = DsDual(0) # current best feasible dual
    primal_f::DsPrimal = DsPrimal() # current best feasible primal
    reopt_candidate_vec::Vector{DsPrimal} = []
    path_hash_set::Set{UInt} = Set{UInt}() # used for storing the path the have been optimized

    # The parameters for the method
    uk::Float64 = 1.0 # The weighting parameter for the quadratic term
    dual_center::DsDual = DsDual(0) # The center of the dual variable 
    kappa::Float64 = 0.5
    J::Vector{Int} = Int[] # The index set
    dual_vec::Vector{DsDual} = DsDual[]
    D_vec::Vector{Float64} = Float64[] # The value of D(λ) for previous iterations
    g_vec = [] # The vector of all previous subgradients
    D_center::Float64 = -Inf
    delta::Float64 = Inf # To measure the progress of the method.
end

