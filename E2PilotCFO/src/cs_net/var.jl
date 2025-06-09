
"""
Decisions on a CsEdge.
"""
@with_kw mutable struct CsTVar{T <: Real}
    tr::T = Inf
    tw::T = Inf
    tc::T = Inf
end

function Base.:(==)(t1::CsTVar, t2::CsTVar)
    return (t1.tr==t2.tr) && (t1.tw==t2.tw) && (t1.tc==t2.tc)
end

tc(t_vec::Vector{T}) where T <: CsTVar = sum([t.tc for t in t_vec])
tw(t_vec::Vector{T}) where T <: CsTVar = sum([t.tw for t in t_vec])
tr(t_vec::Vector{T}) where T <: CsTVar = sum([t.tr for t in t_vec])

function Base.sum(t_vec::Vector{T}) where T <: CsTVar
    t = t_vec[1]
    for i in 2:length(t_vec)
        t+=t_vec[i]
    end
    return t.tr + t.tw + t.tc
end

function Base.:+(t1::CsTVar,t2::CsTVar)
    t1 = deepcopy(t1)
    t1.tr += t2.tr 
    t1.tw += t2.tw
    t1.tc += t2.tc
    return t1
end