"""
PAth selection and Speed Optimization algorithm
Refer to:
Energy-Efficient Timely Transportation of Long-Haul Heavy-Duty Trucks
"""



"""
The PASO algorithm, given source, destination, deadline and vehicle,
return the path & time profile 
output_λ: output λ for performance diagonosis
"""
function paso(net::Network, src::Int64, des::Int64, T::Real; veh::AbstractVehicle, objtype::AbstractObjective = ObjEnergy(), visitonce::Bool=false, breakearly=false, output_step::Bool=true, debug_msg::Bool = true, minimize_t::Bool = false, no_charge::Bool = false, B::Float64 = Inf, xrtol::Float64=1e-3, kwargs...)
    function Delta(lambda::Real)
        
        getdist(net, u, v, lambda) = paso_getdist(net, u, v, lambda, veh, objtype, no_charge, false)
        path0::Vector{Int} = shortest_path(Astar(), net, src, des, lambda=lambda, getdist=getdist, heuristic=heuristic, visitonce=visitonce, breakearly=breakearly) 
        time_v0::Vector{Float64} = [augmentedtime(net, path0[i], path0[i+1], lambda, veh)[1] for i in 1:length(path0)-1]
        total_time = sum(time_v0)
        cost_v::Vector{Float64} = [paso_getdist(net, path0[i], path0[i+1], lambda, veh, objtype, no_charge, true) for i in 1:length(path0)-1]
        tot_cost0::Float64 = sum(cost_v)
        # @debug "Delta called with λ=$lambda, Delta(λ)=$(total_time/3600), desired_T=$(T/3600)"
        return total_time, path0, time_v0, tot_cost0
    end
    # if isinf(B) && isa(veh, EV)
    #     B = veh.cap 
    # end
    if minimize_t && isa(veh, EV)
        func = function (lambda::Real)
            Del_lam = Delta(lambda)[4]
            fval =   Del_lam - B
            # @debug "func called with λ=$lambda, fval=$fval"
            return fval
        end
        # @debug "func for minimize_t"
    else
        func = function (lambda::Real)
            Del_lam = Delta(lambda)[1]
            fval = T - Del_lam
            # @debug "func called with λ=$lambda, fval=$fval"
            return fval
        end
    end
    
    # heuristic = get_heuristic(net, src, des, veh, objtype; visitonce=visitonce)
    heuristic = x -> 0.0
    max_lam = max_lambda(veh) - 1e-6

    method = Roots.A42()
	lam_opt = find_zero(func, (1e-6, max_lam), method; verbose=false, xrtol=xrtol)
    _, path, time_v, tot_cost = Delta(lam_opt*(1+1e-8))
    if minimize_t
        debug_msg && @debug "paso done with lam_opt: $lam_opt, cost=$(tot_cost), desired_cap=$(veh.cap)"
    else
        debug_msg && @debug "paso done with lam_opt: $lam_opt, T=$(sum(time_v)/3600), desired_T=$(T/3600)"
    end
    if output_step
        step = Step(net, veh, path, time_v, objtype; kwargs...)
        return step
    else
        return path, time_v, tot_cost
    end
end


function paso_getdist(net, u::Int, v::Int, lambda::Real, veh::AbstractVehicle, objtype::AbstractObjective, no_charge::Bool, cost_flag::Bool)

    # @timeit g_to "augmenttime" 
    t::Float64, opt_spd::Float64 = augmentedtime(net, u, v, lambda, veh)
    # if isa(veh, EV) && no_charge && is_charge_edge(u, v, net.cs_flag_vec)
    #     return 0.0
    # end
    c = cost(objtype, veh, net, u, v, opt_spd)
    if is_undef_pow(c)
        max_spd = get_max_speed(net, u, v, veh)
        @warn "Got undef cost $c for ($u,$v), $opt_spd, max_spd=$max_spd"
    end
   
    return ifelse(cost_flag, c, c+lambda*t)
end

"""
Create a step from path and time_v
"""
function Step(net::Network, veh::AbstractVehicle, path::Vector{Int64}, time_v::Vector{Float64}, objtype::AbstractObjective = ObjEnergy())
    step = Step()
    node_v = getnode(net, path)
    n_node = length(node_v)
    dis_v = distance3d.(node_v[1:end-1], node_v[2:end])
    div(x,y) = y == 0 ? 0.0 : x/y
    # spd_v = dis_v ./ time_v
    spd_v = div.(dis_v, time_v)
    speedplan = [
        SpeedPlanPoint(node_v[i], spd_v[i], dis_v[i]) for i in 1:n_node-1
    ]
    step.speedplan = speedplan
    cost = total_cost(veh, node_v, spd_v)
    step.summary.cost = cost
    step.summary.duration = sum(time_v)
    step.summary.distance = sum(dis_v)
    return step
end
