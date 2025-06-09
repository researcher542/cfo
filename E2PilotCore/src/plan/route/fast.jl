
"""
Fastest path with minimum time as weight
slow_flag: use the slowest speed for minimum energy consumption
"""
function fastest_path(net::AbsNet, src::Int64, des::Int64; veh::AbstractVehicle=Vehicle(), objtype::AbstractObjective = ObjEnergy(), output_path::Bool = false, slow_flag::Bool = false,  kwargs...)
    function getdist(net, u, v, lambda)
        # dis = distance3d(net, u, v)
        min_t, max_t = minmax_t(net, u, v, veh)
        # spd = get_max_speed(net, u, v, veh)
        # if dis == 0 || spd == 0
        #     return 0.0
        # end
        # t = dis/spd
        return slow_flag ? max_t : min_t
    end
    path = shortest_path(Dijkstra(), net, src, des; getdist=getdist)
    time_v = [getdist(net, path[i], path[i+1], 0) for i in 1:length(path)-1]
    if output_path
        return path, time_v
    end
    set_objtype!(veh, objtype)
    step0 = Step(net, veh, path, time_v, objtype; kwargs...)
    
    return step0
end