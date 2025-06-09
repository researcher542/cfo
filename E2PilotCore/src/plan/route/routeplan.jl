"""
"""



abstract type AbstractShortestPath end
struct Dijkstra <: AbstractShortestPath end
struct Astar <: AbstractShortestPath end
struct SPFA <: AbstractShortestPath end

include("fast.jl")
include("astar.jl")

"""
When we set heurstic to 0 for Astar, it is simply an dijkstra
"""
shortest_path(type::Dijkstra, net::AbsNet, src::Int64, des::Int64 ; lambda::Real = 0.0, getdist::Function=augmentedtime) = shortest_path(Astar(), net, src, des; lambda = lambda, getdist=getdist, heuristic=x->0.0, visitonce=true, breakearly=true)


"""
fastflag: if true, we restrict the speed to maximum speed
"""
function minmax_t(net::Network, i::Int, j::Int, veh::AbstractVehicle, fastflag::Bool = false)
    (i == j) && return (0.0, 0.0)
    dis = distance(net, i, j)
    #max_spd = get_max_speed(net, i, j, ev)
    min_spd,max_spd = get_minmax_speed(net, i, j, veh)
    min_t, max_t = dis/max_spd, dis/min_spd
    if fastflag
        return (min_t, min_t)
    end
    #max_t = max(min_t, dis/g_min_highway_speed) + 1e-3
    return min_t, max_t
end

