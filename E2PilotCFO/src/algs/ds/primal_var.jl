
"""
A subsolution of ds contains a subpath and corresponding scheduling decision

Note that here the charge decision happens at the beginning of the subpath.
"""
@with_kw mutable struct DsSubsolution
    path::Vector{Int} = []
    t_vec::Vector{Float64} = []
    tw::Float64 = 0.0
    tc::Float64 = 0.0
    τ::Float64 = 0.0 # The moment of entering the charging station
    β::Float64 = 0.0 # The init SoC of entering the charging station.
    ics::Int = -1 # The index of charging station in the cs_vec
end

@with_kw_noshow mutable struct DsPrimal
    sub_sol_vec::Vector{DsSubsolution} = []
    beta_cs_vec::Vector{Float64} = []
    tau_cs_vec::Vector{Float64} = []
end

function isvalid(p::DsPrimal)
    if isempty(p.sub_sol_vec)
        return false
    elseif isempty(p.sub_sol_vec[1].t_vec)
        return false
    end
    return true
end

Base.isempty(p::DsPrimal) = isempty(p.sub_sol_vec)

function Base.getindex(primal::DsPrimal, i::Int)
    return primal.sub_sol_vec[i]
end

function Base.show(io::IO, primal::DsPrimal)
    Base.print(io, typeof(primal))
    Base.print(io, "with $(length(primal.sub_sol_vec)) sub_sol_vecs\n")
    for (isol, sol) in enumerate(primal.sub_sol_vec)
        # Base.print(io, "sub_sol_vec[$(isol)]: \n")
        # Base.show(io, sol)
        # Base.print(io, "\n")
    end
    # @show primal.beta_cs_vec 
    # @show primal.tau_cs_vec 
end

function get_travel_time(primal::DsPrimal)
    return sum([sum(sol.t_vec) for sol in primal.sub_sol_vec])
end

function get_wait_time(primal::DsPrimal)
    return sum([sol.tw for sol in primal.sub_sol_vec])
end

function get_charge_time(primal::DsPrimal)
    return sum([sol.tc for sol in primal.sub_sol_vec])
end

get_N(p::DsPrimal) = length(p.sub_sol_vec) - 1
function get_cs_path(primal::DsPrimal)
    path1::Vector{Int} = [
        sol.path[1] for sol in primal.sub_sol_vec
    ]
    path_end::Int = primal.sub_sol_vec[end].path[end]
    push!(path1, path_end)
    # cs_path::Vector{Int} = vcat(path1, path_end)
    cs_path = path1
    return cs_path
end

function get_cs_path_unique(primal::DsPrimal)
    cs_path0 = get_cs_path(primal)
    cs_path = unique!(cs_path0)
    return cs_path
end

function Base.copy(p::DsPrimal)
    # return deepcopy(p)
    return DsPrimal([copy(sol) for sol in p.sub_sol_vec], copy(p.beta_cs_vec), copy(p.tau_cs_vec)) 
end

function hash_path(p::DsPrimal)
    return hash([sol.path for sol in p.sub_sol_vec])
end

"""
convert a primal to path and t_vec
"""
function ds_primal2path(primal::DsPrimal)
    path = reduce(vcat, [sol.path[1:end-1] for sol in primal.sub_sol_vec])
    path = vcat(path, primal.sub_sol_vec[end].path[end])
    t_vec = reduce(vcat, [sol.t_vec for sol in primal.sub_sol_vec])
    return path, t_vec 
end

function get_total_distance(net::Network, primal::DsPrimal)
    path, t_vec = ds_primal2path(primal) 
    dis = sum([distance3d(net, path[i], path[i+1]) for i in 1:length(path)-1])
    return dis
end

function get_total_time(primal::DsPrimal)
    path, t_vec = ds_primal2path(primal) 
    tc_vec = [sol.tc for sol in primal.sub_sol_vec]
    tw_vec = [sol.tw for sol in primal.sub_sol_vec]
    return sum(t_vec) + sum(tc_vec) + sum(tw_vec)
end
