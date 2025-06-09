
"""
compute multiple paso with divided lambda range.
"""

"""
A cache for storing the distance between two nodes with a given lambda.
"""
@with_kw mutable struct PasoDistCache
    # From (src, des, lam) to distance
    dict::Dict{Tuple{Int,Int,Float64}, Float64} = Dict()
    lam_range::Vector{Float64} = []    
    no_charge::Bool = true
    all_pos_lam_vec::Vector{Float64} = [] # The lam inside this vector is guanranteed to be produce postive edges such that we can use Dijkstra to compute the path.
end

"""
For every edge, compute the distance with a given lambda.
"""
function get_paso_dist_cache(net::AbsNet, veh::AbstractVehicle, objtype::ObjEnergy, nlam::Int)
    max_lam = max_lambda(veh) - 1e-6
    lam_range = [x^3 for x in LinRange(0.0, 1.0, nlam)] * max_lam
    dict_pair_vec = Pair{Tuple{Int,Int,Float64}, Float64}[]
    no_charge = true
    
    for lam in lam_range
        for e in Graphs.edges(net.graph)
            (u, v) = Graphs.src(e), Graphs.dst(e)
            dis = paso_getdist(net, u, v, lam, veh, objtype, no_charge, false) 
            push!(dict_pair_vec, (u, v, lam) => dis)
        end
    end
    d = Dict(dict_pair_vec)
    cache = PasoDistCache(;dict=d, lam_range=lam_range, no_charge=no_charge)
    for lam in lam_range
        all_pos_flag = true
        for e in Graphs.edges(net.graph)
            (u, v) = Graphs.src(e), Graphs.dst(e)
            dist = d[(u, v, lam)]
            if dist < 0
                all_pos_flag = false
                # @show u, v, lam, dist
                break
            end
        end
        if all_pos_flag
            push!(cache.all_pos_lam_vec, lam)
        end
    end
    return cache
end

Base.isempty(cache::PasoDistCache) = isempty(lam_range(cache))

function Base.getindex(cache::PasoDistCache, u::Int, v::Int, lam::Float64)
    return cache.dict[(u, v, lam)]
end

"""
"""
function multi_paso(net::Network, src::Int64, des::Int64, veh::AbstractVehicle, objtype::ObjEnergy, n_sol::Int, cache = nothing)
    function Delta(lambda::Real)
        no_charge = true
        function getdist(net, u, v, lambda) 
            if isnothing(cache)
                return paso_getdist(net, u, v, lambda, veh, objtype, no_charge, false)
            else
                return cache[u, v, lambda]
            end
        end
        path0::Vector{Int} = shortest_path(Astar(), net, src, des, lambda=lambda, getdist=getdist, heuristic=heuristic, visitonce=false, breakearly=false) 
        time_v0::Vector{Float64} = [augmentedtime(net, path0[i], path0[i+1], lambda, veh)[1] for i in 1:length(path0)-1]
        total_time = sum(time_v0)
        cost_v::Vector{Float64} = [paso_getdist(net, path0[i], path0[i+1], lambda, veh, objtype, no_charge, true) for i in 1:length(path0)-1]
        tot_cost0::Float64 = sum(cost_v)
        # @debug "Delta called with λ=$lambda, Delta(λ)=$(total_time/3600), desired_T=$(T/3600)"
        return total_time, path0, time_v0, tot_cost0
    end
    
    heuristic = x -> 0.0

    # heuristic = get_heuristic(net, src, des, veh, objtype; visitonce=visitonce)
    if isnothing(cache)
        max_lam = max_lambda(veh) - 1e-6
        # lam_range = LinRange(0.0, max_lam, n_sol)
        _, path, time_v1, tot_cost = Delta(max_lam)
        _, path, time_v2, tot_cost = Delta(0.0)
        min_t = sum(time_v1)
        max_t = sum(time_v2)

        mid_t = min_t + 0.01 * (max_t - min_t)
        # find a lambda such that the time t start to decrease. To avoid redundent computation
        max_lam0 = max_lam
        tot_t = max_t
        while tot_t > mid_t
            max_lam0 = max_lam0 / 2.0
            _, path, time_v, _ = Delta(max_lam0)
            tot_t = sum(time_v)
            if tot_t < mid_t
                # max_lam0 = max_lam0 * 2.0
                break
            end
        end
        lam_range = [x^2 for x in LinRange(0.0, 1.0, n_sol)] * max_lam0
        push!(lam_range, max_lam)
    else
        lam_range = cache.lam_range

        max_lam = maximum(lam_range)
        _, path, time_v1, tot_cost = Delta(max_lam)
        min_t = sum(time_v1)
    end
    # max_lam0 = find_zero(lam -> mid_t - Delta(lam)[1], (0.0, max_lam) ; xrtol=1e-2, verbose=false)
    # @show max_lam0 max_lam
    # @debug "" max_lam0 min_t max_t mid_t
    # @show Delta(max_lam0)[1]

    # @show lam_range
    
    step_vec::Vector{Step} = Step[]

    if !issorted(lam_range)
        lam_range = sort(lam_range)
    end

    t_tol = 1e-3
    for lam in lam_range
        _, path, time_v, tot_cost = Delta(lam)
        step = Step(net, veh, path, time_v, objtype)
        push!(step_vec, step)
        tot_t = sum(time_v)
        if tot_t < min_t + t_tol * abs(min_t)
            break
        end
    end
    lt_func = (s1,s2) -> (s1.summary.duration < s2.summary.duration)
    sort!(step_vec; lt=lt_func, rev=true)


    return step_vec

    # method = Roots.A42()
	# lam_opt = find_zero(func, (0, max_lam), method, verbose=false, xrtol=1e-3)
    # _, path, time_v, tot_cost = Delta(lam_opt*(1+1e-8))
    # if minimize_t
    #     debug_msg && @debug "paso done with lam_opt: $lam_opt, cost=$(tot_cost), desired_cap=$(veh.cap)"
    # else
    #     debug_msg && @debug "paso done with lam_opt: $lam_opt, T=$(sum(time_v)/3600), desired_T=$(T/3600)"
    # end
    # if output_step
    #     step = Step(net, veh, path, time_v, objtype; kwargs...)
    #     return step
    # else
    #     return path, time_v, tot_cost
    # end
end