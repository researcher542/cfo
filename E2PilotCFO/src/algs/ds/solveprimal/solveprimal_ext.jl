

"""
Directly solve the primal problem on stage-expanded graph.
There are N+1 stages
"""
function ds_solve_primal_ext(prob::AbsCfoProb, ext_net::StageExtNetwork, last_primal::DsPrimal, dual::DsDual, option::DsOption, debug::Bool = false)

    (;net, src, des, ev, T, β0) = prob
    (;N, t_mul, b_mul) = option
    cs_vec::Vector{ChargeStation} = prob.net.cs_vec

    @timeit g_to "ds_get_dist_mx" low_dist_mx = ds_get_dist_mx(net, dual, ev, N, option, prob.objtype)
    @timeit g_to "get_cs_cost" cs_cost_vec = ds_get_cs_cost_vec(prob, cs_vec, last_primal, dual, option)

    src = prob.src
    ext_des = orig2ext_idx(ext_net, prob.des, N+1)
    # getdist = let low_dist_mx=low_dist_mx, cs_cost_vec=cs_cost_vec
    #     function getdist(net::StageExtNetwork, u::Int, v::Int, lam) 
    #         return ds_ext_get_dist(net, u, v, low_dist_mx, cs_cost_vec)
    #     end
    # end
    # args = ()
    getdist(net, u, v, lam, low_dist_mx, cs_cost_vec) = ds_ext_get_dist(net, u, v, low_dist_mx, cs_cost_vec); args = (low_dist_mx, cs_cost_vec)

    @timeit g_to "shortest path" path0 = ep.shortest_path(ep.Astar(), ext_net, src, ext_des; 
        getdist=getdist, getdistargs=args, visitonce=true, breakearly=true)
    w_vec0 = [getdist(ext_net, path0[i], path0[i+1], 0.0, args...) for i in 1:length(path0)-1 ]
    path_obj0 = sum(w_vec0)
    path1 = [ ext2orig_idx(ext_net, u) for u in path0] 

    ### Recover the primal variable from the path on the ext graph.
    sub_sol_vec = DsSubsolution[]
    subpath = Int[]
    sub_t_vec = Float64[]

    tw = 0.0
    tc = 0.0
    τ = 0.0
    β = prob.β0
    ics = 0
    sigma_obj = 0.0
    w_obj = 0.0
    for ipath1 in 1:length(path1)-1
        (u, istage_u) = path1[ipath1]
        (v, istage_v) = path1[ipath1+1]
        push!(subpath, u)
        if istage_v == istage_u + 1
            # If this is the next stage.
            @assert (u == v)
            if isempty(sub_t_vec)
                push!(sub_t_vec, 0.0) 
            end
            if length(subpath) == 1
                push!(subpath, subpath[1])
            end
            subsol = DsSubsolution(subpath, sub_t_vec, tw, tc, τ, β, ics)
            push!(sub_sol_vec, subsol)
            subpath = Int[]
            sub_t_vec = Float64[]

            #
            ics = findfirst(cs -> (cs.idx == u), cs_vec)
            cs0 = cs_vec[ics]
            # ta = prob.carbon_dict[cs0.region]
            ta = get_price_ta(prob, cs0.idx, prob.predict_mode)
            cost, x = ds_get_cs_cost(prob, last_primal, dual, istage_u, cs0.idx, ta, option)
            tc, tw, β, τ = x
            # @show istage_u, cost x
            sigma_obj += cost
        else
            # If u v is in the same stage.
            @assert istage_u == istage_v
            (c, t) = ds_get_dist(net, istage_u, u, v, dual, ev, option, prob.objtype)
            w_obj += c
            push!(sub_t_vec, t)
        end
        
    end

    if isempty(subpath)
        ## If the last edge is (des, N) -> (des, N+1)
        subpath = Int[prob.des, prob.des]
        sub_t_vec = Float64[0.0]
    else
        push!(subpath, prob.des)
    end
    subsol = DsSubsolution(subpath, sub_t_vec, tw, tc, τ, β, ics)
    push!(sub_sol_vec, subsol)

    # push!(subpath, prob.des)
    # push!(sub_t_vec, 0.0)
    # subsol = DsSubsolution(subpath, sub_t_vec, tw, tc, τ, β, ics)
    # push!(sub_sol_vec, subsol)

    D1b = β0*dual.μ[1] * b_mul
    DN1t = dual.λ[N+1] * t_mul * T
    DN1b = dual.μ[N+1] * b_mul * option.β_lb

    # @show sigma_obj w_obj

    # @show sigma_obj + w_obj
    # @show sigma_obj w_obj DN1b -D1b -DN1t
    # @show sigma_obj + w_obj + DN1b - D1b - DN1t
    # @show path_obj0 w_vec0

    primal = DsPrimal(;sub_sol_vec=sub_sol_vec)
    return primal
end


