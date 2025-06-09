"""
Simulate the solution to validate.
"""

convert_t_vec(t_vec) = t_vec
function convert_t_vec(t_vec::Vector{Float64})  
    charge_time = 80.0 * 60
    n_node = length(t_vec)
    t_vec_new = [
    CsTVar(
        t_vec[i], 
        i==n_node ? 0.0 : g_min_wait_time, 
        i==n_node ? 0.0 : charge_time
        ) for i in 1:n_node
    ]
    return t_vec_new
end

"""
simulate over cs_net

charge: if we need to charge during the path

charge_first, if true, on a cs_edge, we perform the decision tw,tc,tr.
    If false, we perform the decisions tr,tw,tc
"""
function simulate_cs_net(prob::AbsCfoProb, path, t_vec0; check_flag::Bool=true, charge_first::Bool = false)
    n_node = length(path)
    t_vec = convert_t_vec(t_vec0)
    (;net, src, des, ev, β0, start_time, T) = prob
    @assert net.cs_flag_vec[src]
    charge_time = 0.0 * 60
    next_node = outneighbors(net, src)[1]
    if charge_first
        path_all, t_vec_all = Int[src, next_node], Float64[charge_time]
    else
        path_all, t_vec_all = Int[src], Float64[]
    end
    
    for inode in 1:n_node-1
        i = path[inode]
        j = path[inode+1]
        tvar::CsTVar = t_vec[inode]
        (;tr,tw,tc) = tvar
        # if !net.cs_flag_vec[i]
        if charge_first 
            # This is for the decision order tw,tc,tr
            i1 = outneighbors(net, i)[1]
            i_pre = inneighbors(net, i)[1]
            j_pre = inneighbors(net, j)[1]
            path_paso, t_vec_paso, tot_cost = paso(net, i1, j_pre, tr; veh=ev, output_step=false, no_charge=true, visitonce=true, debug_msg=false)
            path_all = vcat(path_all, i, path_paso)
            t_vec_all = vcat(t_vec_all, tw, tc, t_vec_paso)
        else
            # This is for the decision order tr,tw,tc
            i1 = outneighbors(net, i)[1]
            path_paso, t_vec_paso, tot_cost = paso(net, i1, j, tr; veh=ev, output_step=false, no_charge=true, visitonce=true, debug_msg=false)
            # @debug tot_cost
            next_node = outneighbors(net, j)[1]
            path_all = vcat(path_all, path_paso[2:end], next_node)
            t_vec_all = vcat(t_vec_all, t_vec_paso[1:end-1], tw, tc)
        end
    end

    # ret = (simulate(net, path_all, t_vec_all, ev, β0, start_time, T; kwargs...)..., t_vec)
    return simulate(prob, path_all, t_vec_all; check_flag=check_flag)
end


function simulate(prob::AbsCfoProb, path, t_vec; check_flag::Bool = true)
    (;net, ev, β0, start_time, T) = prob
    return simulate(net, path, t_vec, ev, β0, start_time, T; check_flag=check_flag)
end

# function simulate(net::Network, path::Vector{Int}, t_vec::Vector, ev::EV, β0::Float64, start_time::DateTime, T::Real; check_flag::Bool = true)
#     #(;t,β,τ,βp,τp) = primal
#     obj::Float64 = 0.0
#     # cs_vec = net.cs_vec
#     cs_flag_vec = net.cs_flag_vec
#     n_node = length(path)
#     β_vec = zeros(Float64, n_node)
#     τ_vec = zeros(Float64, n_node)
#     β_vec[1] = β0
#     τ_vec[1] = 0.0
#     infeasible_flag = false
#     for inode in 1:n_node-1
#         τ = τ_vec[inode]
#         β = β_vec[inode]
#         if check_flag && ( (β < -1e-2 || β > ev.cap+1e-2) || (τ<0 || τ > T+1e-2) )
#             @warn "infeasible plan at i=$inode" β τ/3600 T/3600
#             check_flag = false
#             infeasible_flag = true
#         end
#         i = path[inode]
#         j = path[inode+1]
#         #tij = primal.t[i,j]
#         tij = t_vec[inode]
#         @assert has_edge(net.graph, i, j)  
#         road_type = get_road_type(i, j, cs_flag_vec)
#         β_vec[inode+1] = get_next_b(net, i, j, tij, β, ev; no_clamp=false)
#         # β - Δsoc(net, i, j, tij, β, ev, road_type)
#         τ_vec[inode+1] = tij + τ
#         obj += cfo_objective(road_type, β, tij, τ, start_time, ev.cap, price_ta)
#     end
#     infeasible_flag && @warn @sprintf "Infeasible with min_β=%.3f %%" minimum(β_vec)/ev.cap*100
#     return obj, β_vec, τ_vec
# end