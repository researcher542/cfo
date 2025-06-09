"""
An dict that behaves like an array
"""

@with_kw mutable struct ArrayDict{T}
    d::Dict{Tuple{Int, Int}, T} = Dict{Tuple{Int, Int}, Int}()
end

function Base.getindex(d::ArrayDict, i::Int, j::Int)
    return d.d[(i,j)] 
end

function Base.setindex!(d::ArrayDict{T}, val::T, i::Int, j::Int) where T
    d.d[i,j]  = val
end