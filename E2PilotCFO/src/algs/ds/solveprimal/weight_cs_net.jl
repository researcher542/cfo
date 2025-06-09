
"""
Get the distance matrix for each pair of charging stations.
Use charging station graph (CsNet) as the low level distance matrix
This is the one that should be called by others.
mx[i][u,v] is the SP^i(u,v)
"""
function ds_get_cs_sp_mx(prob::AbsCfoProb, cs_vec::Vector{ChargeStation}, dual::DsDual, N::Int, option::DsOption{DsCsNet}) 
    (;net, src, des, ev) = prob
    cs_net = net.cs_net

    n_node = nv(net)

    # append src,des in the node list
    cs_idx_vec = vcat(src, des, [cs.idx for cs in cs_vec])
    I = [(i, j)[1] for i in cs_idx_vec for j in cs_idx_vec]
    J = [(i, j)[2] for i in cs_idx_vec for j in cs_idx_vec]
    cs_dist_mx = [spzeros(n_node, n_node) for _ in 1:N+1]
    # Threads.@thread 
    for i in 1:N+1
        V = zeros(length(I)) 
        fill!(V, Inf)

        Threads.@threads for isrc in eachindex(cs_idx_vec)
            for (ides,des) in enumerate(cs_idx_vec)
                src1 = cs_idx_vec[isrc]
                des1 = cs_idx_vec[ides]
                ivv = (isrc-1) * (length(cs_idx_vec)) + ides
                edge = get_edge(cs_net, src1, des1)
                if isnothing(edge)
                    V[ivv] = Inf
                else
                    cost, t = ds_cs_net_get_dist(edge, dual, i, option)
                    V[ivv] = cost
                end
                # cs_dist_mx[i][src, des] = dist_mx_i[des]
                # @show isrc ides iv ivv
                # @assert iv == ivv 
                # iv += 1
            end
        end
        
        AA = sparse(I, J, V, n_node, n_node)
        cs_dist_mx[i] = AA
    end
    return cs_dist_mx
end


function ds_cs_net_get_dist(edge::CsEdge, dual::DsDual, i::Int, option::DsOption{DsCsNet})
    μi = dual.μ[i]
    λi = dual.λ[i]
    (;min_t, max_t) = edge
    (;t_mul,b_mul) = option
    t_vec = [x[1] for x in edge.t_cost_vec]
    cost_vec = [x[2] for x in edge.t_cost_vec]
   

    if λi == 0 && μi == 0
        # t::Float64 = (min_t+max_t)/2
        t::Float64 = min_t
        cost::Float64 = 0.0
    elseif μi == 0
        @assert(λi > 0)
        t = min_t
        cost = λi*t*t_mul
    else
        if length(t_vec) == 1
            t_vec = [t_vec[1], t_vec[1]]
            cost_vec = [cost_vec[1], cost_vec[1]]
        end
        pwl = PWL(;x=t_vec, y=cost_vec)
        function f_obj(t)
            e_cost = pwl(t, true)   
            cost = e_cost*μi*b_mul + λi*t*t_mul
            return cost
        end
        if (min_t == max_t)
            cost = f_obj(min_t)
            t = min_t
        else
            # @show min_t max_t
            @assert (min_t <= max_t)
            cost, t = golden_section(f_obj, min_t, max_t)
        end
        
        if option.pos_cost && cost < 0
            cost = 0.0
        end
    end

    return (cost, t)
end