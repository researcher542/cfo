

abstract type AbsDualState end

"""
"""
@with_kw mutable struct DsState <: AbsDualState
    k::Int = 1
    r::Int = 0
   
    beta_k::Int = 1
    v::Int = 0
    phaseI::Bool = true

    # debug information
    cons_vio::Float64 = 0.0 # The constraint violation
    obj::Float64 = 0.0 # The raw objective
    lag_obj::Float64 = 0.0 # the Lagrangian objective
    dual::DsDual = DsDual(0)
    primal::DsPrimal = DsPrimal()
    stepsize::Float64 = 0.0

    # The modified sub gradient
    grad_t_m::Vector{Float64} = zeros(0)
    grad_b_m::Vector{Float64} = zeros(0)

    # _c stands for the current best Lagrangian objective
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

    # The unchanged parameters
    theta0::Float64 = 1e20 # The upper bound of the objective
    r1::Int = 5
    eps0::Float64 = 1e-3
    v_ub::Int = 3
    beta_ub::Int = 120 # Note that this is NOT the SoC, but a parameter in the algorithm to update the step size
    eps::Float64 = 1e-3

    state_lock = ReentrantLock()
end


function Base.getproperty(s::AbsDualState, sym::Symbol)
    if sym == :λ
        return s.dual.λ
    elseif sym == :μ
        return s.dual.μ
    elseif sym == :β
        # return the first charging staion for debug...
        return [sol.β/3.6e6 for sol in s.primal.sub_sol_vec]
    elseif sym == :τ
        return [sol.τ/3600.0 for sol in s.primal.sub_sol_vec]
    elseif sym == :tc
        return [sol.tc/3600.0 for sol in s.primal.sub_sol_vec]
    else
        return Base.getfield(s, sym)
    end
end

