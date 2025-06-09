
"""
Directly solve the problem that 
min c(t) + λt
The solution here should be continuous w.r.t. λ
"""
function augmentedtime_raw(net::Network, u::Int, v::Int, lambda::Real, veh::AbstractVehicle)
    # @debug "augmenttime_raw called."
    opt_spd = optimal_speed_raw(net, u, v, lambda, veh) 
    dis = distance(net, u, v)
    t = dis / opt_spd
    return t, opt_spd
end

"""
c(t) = tf(d/t) + λ t
dc(t)/dt = f(v) - v*f(v)
We directly solve the problem min c(t) with golden_section method
"""
function optimal_speed_raw(net::Network, u::Int, v::Int, lambda::Real, veh::AbstractVehicle)
    min_spd, max_spd = get_minmax_speed(net, u, v, veh)
    theta = grade(net, u, v)
    dis = distance(net, u, v)
    function f_obj(spd)
        t = dis / spd
        fval = output_power(veh, theta, spd)
        obj = lambda * t + fval*t
        return obj
    end
    
    _, opt_spd = golden_section(f_obj, min_spd, max_spd; show_trace=false)

    return opt_spd 
end