
@with_kw mutable struct StageExtNetwork{G} <: AbsNet
    g::G = DiGraph()
    N::Int = 0 # Number of charging stops, N+1 is the number of stages.
    orig_nv::Int = 0
    orig_ne::Int = 0
    cs_idx_vec::Vector{Int} = zeros(0)
    orig_des::Int = -1
    orig_src::Int = -1
    out_neigh_dict::OrderedDict{Int, Vector{Int}} = OrderedDict{Int, Vector{Int}}() # The dict to store the buffer of outneighbors in the stage ext network.
end

Graphs.nv(net::StageExtNetwork) = (net.N+2) * Graphs.nv(net.g)
Graphs.ne(net::StageExtNetwork) = (net.N+2) * Graphs.ne(net.g)

function orig2ext_idx(net::StageExtNetwork, u_orig::Int, istage::Int)::Int
    return (istage - 1) * net.orig_nv + u_orig
end

function ext2orig_idx(net::StageExtNetwork, u_ext::Int)
    nv = net.orig_nv
    istage = u_ext ÷ nv + 1
    u_orig = mod(u_ext, nv)
    if u_orig == 0
        istage -= 1
        u_orig = nv
    end
    return (u_orig, istage)
end

function Graphs.outneighbors(net::StageExtNetwork, u_ext::Int)
    if haskey(net.out_neigh_dict, u_ext)
        return net.out_neigh_dict[u_ext]
    end

    (u_orig, istage) = ext2orig_idx(net, u_ext)
    if (istage == net.N + 2) 
        return Int[]
    end

    neigh_vec::Vector{Int} = Graphs.outneighbors(net.g, u_orig)
    # neigh_vec1::Vector{Int} = [orig2ext_idx(net, nei, istage) for nei in neigh_vec]
    # neigh_vec1::Vector{Int} = Int[]
    # @inbounds for nei::Int in neigh_vec
    #     push!(neigh_vec1, orig2ext_idx(net, nei, istage))
    # end
    neigh_vec1::Vector{Int} = map(nei -> orig2ext_idx(net, nei, istage), neigh_vec) 
    if ((u_orig in net.cs_idx_vec) && (istage <= net.N)) || (u_orig == net.orig_des)
        ## If u is the charging station or the destination, add an virtual edge
        next_u  = orig2ext_idx(net, u_orig, istage + 1)
        if u_orig != net.orig_src
            ## We ignore the case where the src is a charging station.
            push!(neigh_vec1, next_u) 
        end
    end

    net.out_neigh_dict[u_ext] = neigh_vec1

    return neigh_vec1
end

function ds_construct_stage_ext_graph(prob::AbsCfoProb, option::DsOption)
    (;net, src, des) = prob
    (;N) = option
    cs_vec::Vector{ChargeStation} = prob.net.cs_vec
    cs_idx_vec = [cs.idx for cs in cs_vec]
    orig_nv = nv(net.g)
    orig_ne = ne(net.g)
    
    return StageExtNetwork(;
        g=deepcopy(net.g), N=N, orig_ne=orig_ne, orig_nv=orig_nv, cs_idx_vec=cs_idx_vec, orig_des = des, orig_src=src,
    )
end

function ds_ext_get_dist(net::StageExtNetwork, u::Int, v::Int, low_dist_mx, cs_cost_vec)::Float64
    (u_orig, u_istage) = ext2orig_idx(net, u)
    (v_orig, v_istage) = ext2orig_idx(net, v)
    if (u_orig == v_orig)
        @assert (u_istage + 1 == v_istage)
        if u_orig == net.orig_des
            return 0.0
        elseif u_orig in net.cs_idx_vec
            ics = findfirst(idx-> idx == u_orig, net.cs_idx_vec)
            return cs_cost_vec[u_istage][ics]
        end
    else
        @assert (u_istage == v_istage) 
        # @debug "" u u_orig u_istage v v_orig v_istage
        dist = low_dist_mx[u_istage][u_orig, v_orig]
        return dist
    end
    @assert false "should not reach here."
    return NaN
end