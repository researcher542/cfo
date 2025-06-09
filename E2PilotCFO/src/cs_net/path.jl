
@with_kw mutable struct CsPathState
    min_cost::Float64 = 0.0
    max_cost::Float64 = 0.0
    min_t::Float64 = 0.0
    max_t::Float64 = 0.0
    dis::Float64 = 0.0
end

function Base.:+(s1::CsPathState, s2::T) where T <: Union{CsPathState, CsEdge}
    s3 = CsPathState(;
        min_cost = s1.min_cost + s2.min_cost,
        max_cost = s1.max_cost + s2.max_cost,
        min_t = s1.min_t + s2.min_t,
        max_t = s1.max_t + s2.max_t,
        dis = s1.dis + s2.dis
    )
    # s3 = deepcopy(s1)
    # s3.min_cost += s2.min_cost 
    # s3.max_cost += s2.min_cost 
    # s3.min_t += s2.min_t 
    # s3.max_t += s2.max_t
    # s3.dis += s2.dis
    return s3
end

function Base.:-(s1::CsPathState, s2::T) where T <: Union{CsPathState, CsEdge}
    s3 = CsPathState(;
        min_cost = s1.min_cost - s2.min_cost,
        max_cost = s1.max_cost - s2.max_cost,
        min_t = s1.min_t - s2.min_t,
        max_t = s1.max_t - s2.max_t,
        dis = s1.dis - s2.dis
    )
    # s3 = deepcopy(s1)
    # s3.min_cost -= s2.min_cost 
    # s3.max_cost -= s2.min_cost 
    # s3.min_t -= s2.min_t 
    # s3.max_t -= s2.max_t 
    # s3.dis -= s2.dis
    return s3
end

@with_kw mutable struct CsPath
    path::Vector{CsEdge{Float64}} = []
    path_int::Vector{Int} = []
    state::CsPathState = CsPathState()
end

Base.length(p::CsPath) = length(p.path)

function push(path::CsPath, edge::CsEdge)
    path0 = deepcopy(path)
    return push!(path0, edge)
end

function Base.push!(path::CsPath, edge::CsEdge)
    if isempty(path.path_int)
        push!(path.path_int, edge.src)
    end
    @assert path.path_int[end] == edge.src
    push!(path.path, edge) 
    push!(path.path_int, edge.des)
    path.state = path.state + edge 
    return path
end

function Base.pop!(path::CsPath)
    edge = path.path[end] 
    path.state -= edge
    pop!(path.path)
    pop!(path.path_int)
end



