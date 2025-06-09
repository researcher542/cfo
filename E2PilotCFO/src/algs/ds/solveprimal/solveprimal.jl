include("weight.jl")
include("weight_cs_net.jl")
include("stage_ext_net.jl")
include("solveprimal_ext.jl")


function ds_solve_primal(prob::AbsCfoProb, last_primal::DsPrimal, dual::DsDual, option::DsOption, debug::Bool = false)
    (;net, src, des, ev, T, β0) = prob
    (;N, t_mul, b_mul) = option
    cs_net = net.cs_net
    cs_vec::Vector{ChargeStation} = prob.net.cs_vec
    # low_dist_mx = ds_get_dist_mx(net, dual, ev, N)
    @timeit g_to "get_cs_sp_mx" cs_sp_mx = ds_get_cs_sp_mx(prob, net.cs_vec, dual, N, option)
    # @debug "cs_sp_mx got." # [Matrix(x) for x in cs_sp_mx]
    @timeit g_to "get_cs_cost" cs_cost_vec = ds_get_cs_cost_vec(prob, cs_vec, last_primal, dual, option)
    # @debug "cs_cost_vec got." # cs_cost_vec

    n_cs = length(cs_vec)
    g_score = fill(Inf, N, n_cs)
    parents = ArrayDict()
	open_set = HeapPriorityQueue{Tuple{Int, Int}, Float64}()
    enqueue!(open_set, (0, 0), 0.0)
    des_score = [Inf]
    function get_score(i::Int, ics::Int)
        if i == 0 
            return 0.0
        elseif i == N+1
            return des_score[1]
        else
            return g_score[i, ics]
        end
    end
    function set_score!(i::Int, ics::Int, s::Float64)
        @assert i != 0
        if i == N+1
            des_score[1] = s
        else
            g_score[i, ics] = s
        end
    end

    while !isempty(open_set)
        (i,ics) = dequeue!(open_set)

        # The destination is reached
        if i == N+1
            break
        end
        cs_cost = (i == 0) ? 0.0 : cs_cost_vec[i][ics]
        u_score = get_score(i, ics)
        for ics_next in eachindex(cs_vec)
            # (ics == ics_next) && continue
            (cs_vec[ics_next].idx == src) && continue
            # (cs_vec[ics_next].idx == des) && continue
            u = (i == 0) ? src : cs_vec[ics].idx
            v = cs_vec[ics_next].idx

            if (i == N)
                v = des
            end

            # cs_edge = get_edge(cs_net, u, v)
            if u == v
                dist = 0.0
                cs_cost = 0.0
                weight = 0.0
            elseif !has_edge(cs_net, u, v)
                dist = Inf
                continue
            else
                weight::Float64 = cs_sp_mx[i+1][u,v]
                dist::Float64 = cs_cost + weight
            end
            tentative_g_score = u_score + dist

            if tentative_g_score < get_score(i+1, ics_next)
                set_score!(i+1, ics_next, tentative_g_score)
                # @debug "updating score" (i+1) ics ics_next tentative_g_score u_score dist
                parents[i+1, v] = u
                enqueue!(open_set, (i+1, ics_next), tentative_g_score)
            end

        end
    end

    # @show get_score(N+1, des)
    # recover the path with recharging stops
    cs_path = zeros(Int, N+2)
    cs_path[N+2] = des
    # @show prob.src prob.des prob.T prob.β0 prob.objtype
    if !haskey(parents.d, (N+1, des))
        @show prob.src, prob.des
    end
    ## If no path here, maybe because the problem is infeasible.
    ## We preprocess the netowrk to get the lb,ub of tau for each charging stations.
    cs_path[N+1] = parents[N+1, des]
    for i in N+1:-1:1
        cs_path[i] = parents[i, cs_path[i+1]]
    end
    @assert cs_path[1] == src
    # @show cs_path

    for i in 1:N
        @assert (has_edge(cs_net, cs_path[i], cs_path[i+1]) || cs_path[i] == cs_path[i+1]) 
    end

    # from cs_path to decisions
    sub_sol_vec = DsSubsolution[]
    sigma_obj = 0.0
    w_obj = 0.0
    for i in 1:N+1

        if cs_path[i] == cs_path[i+1]
            path = [cs_path[i], cs_path[i+1]]
            t_vec = [0.0]
        else
            path, t_vec, w_cost = ds_recover_primal_edge(prob, cs_path, i, dual, option) 
            w_obj += w_cost
        end

        if i == 1
            β = prob.β0
            τ,tc,tw = 0.0, 0.0, 0.0
            ics = 0
        else
            ics_vec = findall(cs->(cs.idx == cs_path[i]), cs_vec)
            @assert(length(ics_vec) == 1)
            ics = ics_vec[1]
            cs0 = cs_vec[ics]
            ta = net.carbon_dict[cs0.region]
            cost, x = ds_get_cs_cost(prob, last_primal, dual, i-1, cs0.idx, ta, option)
            tc, tw, β, τ = x
            sigma_obj += cost
        end
        subsol = DsSubsolution(path, t_vec, tw, tc, τ, β, ics)
        push!(sub_sol_vec, subsol)
        
    end

    D1b = β0*dual.μ[1] * b_mul
    DN1t = dual.λ[N+1] * t_mul * T
    DN1b = dual.μ[N+1] * b_mul * option.β_lb
    # @show "before" sigma_obj + w_obj
    @assert (sigma_obj <= 1e19) (@exfiltrate; @show prob; "No feasible sigma?")
    @show sigma_obj w_obj DN1b -D1b -DN1t
    @show sigma_obj + w_obj + DN1b - D1b - DN1t

    return DsPrimal(;sub_sol_vec=sub_sol_vec)
end

"""
recover the primal solution on a sub-path. The underlying graph is the original graph.
"""
function ds_recover_primal_edge(prob::AbsCfoProb, cs_path::Vector{Int}, i::Int,  dual::DsDual, option::DsOption{DsNetwork} )
    (;net, ev) = prob
    getdist = (net, u, v, lam) -> ds_get_dist(net, i, u, v, dual, ev, option, prob.objtype)
    path = shortest_path(Astar(), net, cs_path[i], cs_path[i+1]; getdist=getdist)
    cost_t_vec = [getdist(net, path[ipath], path[ipath+1], 0.0) for ipath in 1:length(path)-1]
    t_vec = [x[2] for x in cost_t_vec]
    cost_vec = [x[1] for x in cost_t_vec]
    w_cost = sum(cost_vec)
    return path, t_vec, w_cost 
end

function ds_recover_primal_edge(prob::AbsCfoProb, cs_path::Vector{Int}, i::Int,  dual::DsDual, option::DsOption{DsCsNet})
    (;net, ev) = prob
    cs_net = net.cs_net
    src1 = cs_path[i]
    des1 = cs_path[i+1]
    μi = dual.μ[i]
    λi = dual.λ[i]
    (;t_mul,b_mul) = option
    edge = get_edge(cs_net, src1, des1)
    if isnothing(edge)
        T = Inf
    else
        (w_cost0, T) = ds_cs_net_get_dist(edge, dual, i, option)
    end
    (path, t_vec, e_cost) = paso(net, src1, des1, T; veh=ev, debug_msg = false, output_step = false)
    tot_t = sum(t_vec)
    w_cost = e_cost*μi*b_mul + λi*tot_t*t_mul
    return path, t_vec, w_cost 
end