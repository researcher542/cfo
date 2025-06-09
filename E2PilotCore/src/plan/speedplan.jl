

"""
Given lambda, caculate the corresponding required time

return: required time 
        corresponding decision [ speed ]
"""
function augmentedtime(veh::AbstractVehicle, lambda::Real, theta::Float64, dis::Float64, max_spd::Float64)
    opt_spd = optimal_speed(veh.sol_grid, theta, lambda)
    # max_spd = request_max_spd(node1)
    opt_spd = min(opt_spd, max_spd)
    opt_time = dis/opt_spd
	return (opt_time, opt_spd)
end

#augmentedtime(net::Network, way::Way, lambda::Real, veh::Vehicle) = augmentedtime(net, way, lambda, veh.sol_grid)

"""
"""
function augmentedtime(net::Network, way::AbstractWay, lambda::Real, veh::AbstractVehicle)
    # Optimized for the case when there is only two nodes
    max_spd = get_max_speed(net, way.src, way.des, veh)
    sol_grid = veh.sol_grid
    # max_spd = g_max_highway_speed
    if length(way.nodes) == 2
        # node1 = getnode(net, way.src)
        # node2 = getnode(net, way.des)
        theta = way.grades[1]
        dis = way.distances[1]
        return augmentedtime(veh, lambda, theta, dis, max_spd)
	    #return augmentedtime(node1, node2, lambda, sol_grid; grade=theta, dis, max_spd=max_spd)
    end
	nodes = getnode(net, way.nodes)
	return augmentedtime(nodes, lambda, sol_grid; grade_v=way.grades, dis_v=way.distances, max_spd=max_spd)
end

augmentedtime(net::Network, src::Int64, des::Int64, lambda::Real, veh::AbstractVehicle) = augmentedtime(net, getway(net, src, des), lambda, veh)


