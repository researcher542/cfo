
@with_kw mutable struct DsDual
    λ::Vector{Float64} # For the time related constraint
    μ::Vector{Float64} # For the SoC lower bound constraint
    # μu::Vector{Float64} # For the SoC upper bound constraint
end


"""
N: number of recharging stops.
"""
function DsDual(N::Int, default::Float64 = 1e-3)
    # return DsDual(zeros(N+1), zeros(N+1), zeros(N)) 
    dual = DsDual(zeros(N+1), zeros(N+1)) 
    fill!(dual.λ, default)
    fill!(dual.μ, default)
    return dual
end

function ds_dual2vec(dual::DsDual)
    return vcat(dual.λ, dual.μ)
end

function ds_vec2dual(vec)
    Np1 = length(vec) ÷ 2
    dual = DsDual(vec[1:Np1], vec[Np1+1:end])
    return dual
end