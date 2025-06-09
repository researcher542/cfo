
@with_kw_noshow mutable struct DsBnBBox{T <: AbstractArray}
    xlb::T #tc,tw,b,tau
    xub::T
    flb::Float64 = -Inf
    fub::Float64 = Inf
    div_idx::Int = 1
end

function Base.show(io::IO, box::DsBnBBox{T}) where T
    Base.print(io, typeof(box))
    Base.print(io, "\n")
    for sym in [:xlb, :xub, :flb, :fub, :div_idx]
        Base.print(io, "    " *string(sym) * ": ")
        Base.show(io, getproperty(box, sym))
        Base.print(io, "\n")
    end
    # Base.print(io, "Solution grid with spd_range $(sg.spd_range).\n theta_range $(sg.theta_range)")
    # Base.show(io, box.xlb)
    # Base.show(io, box.xub)
end



function volume(box::DsBnBBox)
    x_diff = (box.xub-box.xlb)
    vol = reduce(*, x_diff) 
    return vol 
end

function Base.copy(box::DsBnBBox)
    return DsBnBBox(Base.copy(box.xlb), Base.copy(box.xub), box.flb, box.fub, box.div_idx) 
end

function Base.lt(o::Base.Order.ForwardOrdering, b1::DsBnBBox{T}, b2::DsBnBBox{T}) where T
    return  (b1.flb < b2.flb)  ||
            ( (b1.flb == b2.flb) && b1.xlb[1] < b2.xlb[1]) 
            # ( (b1.flb == b2.flb) && b1.xlb[1] == b2.xlb[1] && b1.div_idx < b2.div_idx)
end

"""
A branch and bound method to solve compute the value of σ, i.e., the cost of a charging station.
istage: the i-th charging stops (stage).


return (cost, x)
"""
function ds_get_cs_cost_bnb(prob::AbsCfoProb, primal::DsPrimal, dual::DsDual, istage::Int, cs_road_idx::Int, ta::TS.TimeArray, option::DsOption)

    B = prob.ev.cap
    st = prob.start_time
    # tc, tw, b, tau
    lb, ub = ds_get_cs_var_range(prob, primal, istage, cs_road_idx, option)
    if ub[4] < lb[4]
        # cost = ds_sigma_objective(lb, primal, dual, istage, ta, st, B, option)
        cost = ds_sigma_objective(prob, lb, primal, dual, istage, cs_road_idx, option)
        return (cost + 1e20, lb)
    end
    ta1 = TS.from(ta, st - Hour(10))
    ta2 = TS.to(ta1, st + τ2second(prob.T) + Hour(10))
    # From box to lower bound
    pq = HeapPriorityQueue{DsBnBBox{MVector{4, Float64}}, DsBnBBox{MVector{4, Float64}}}()
    # pq = PriorityQueue{DsBnBBox{MVector{4, Float64}}, Float64}()

    box =  DsBnBBox(;xlb=lb, xub=ub, div_idx=1)
    box = ds_bnb_update_lbub!(box, primal, dual, istage, ta2, st, B, option)

    enqueue!(pq, box, box)
    fub = box.fub
    tol = 5e-3
    best_box = box
    cur_box = box
    # for iter in 1:10000
    rel_diff = Inf
    abs_diff = Inf
    # idx_cnt_vec = zeros(1000)
    remain_vol = volume(box)
    tot_vol = volume(box)
    iter::Int = 0
    while iter <= 1_000_000
        iter += 1
    # for iter in 1:10
        isempty(pq) && break
        cur_box = dequeue!(pq)
        # @show iter cur_box
        if iter % 1_000_00_000 == 0
        # if iter % 1 == 0
            @show iter cur_box best_box
            # @show cur_box.xlb cur_box.xub cur_box.fub cur_box.flb cur_box.div_idx
            # @show best_box.xlb best_box.xub best_box.fub best_box.flb 
            # @show best_box.div_idx
            # @show cur_box.flb
            @show fub diff
            @show length(pq) diff
            @show remain_vol/tot_vol
            println()
        end
        # println()
        abs_diff = (fub-cur_box.flb)
        rel_diff = abs_diff/max(abs(fub), abs(cur_box.flb))
        if rel_diff  < tol || abs_diff < tol
            # @show diff iter remain_vol/tot_vol
            break
        end
        
        # divide box
        box1, box2 = divide(cur_box)
        # update the divided box and global upper bound
        for box_i in SA[box1, box2]
            box_i = refine!(box_i, B)
            box_i = ds_bnb_update_lbub!(box_i, primal, dual, istage, ta2, st, B, option)
            if box_i.fub < fub
                fub = box_i.fub
                best_box = copy(box_i)
            end
            if box_i.flb == cur_box.flb
                # idx_cnt_vec[box_i.div_idx] += 1
                # @show box_i.div_idx
            end
            @assert box_i.flb + 1e-12  >= cur_box.flb (@exfiltrate; "divide should provide larger lower bound")
            # enqueue the box only when it is not be pruned
            # enqueue!(pq, box_i, box_i) 
            if (box_i.flb <= fub) # && (box_i.flb != box_i.fub)
                enqueue!(pq, box_i, box_i) 
            elseif box_i.flb ≈ fub
                # diff = 0.0
                enqueue!(pq, box_i, box_i) 
                # break
            else
                remain_vol -= volume(box_i)
            end
        end
    end
    @assert (abs_diff <= 5e-2 || rel_diff <= 5e-2) (@exfiltrate; "diff larger than tolerance")
    # x = (best_box.xlb+best_box.xub)/2
    # This is biased, xlb seems likely to achieve lower bound.
    # x, val = best_box.xlb
    # @show diff
    x, val = ds_bnb_get_ub(best_box, primal, dual, istage, ta, st, B, option)
    cost = best_box.fub
    # @show iter
    # @show diff remain_vol
    # @show best_box.xlb best_box.xub fub

    # @debug "" cs_road_idx cost best_box primal dual
    # @show lb ub best_box
    return (cost, x)
end

"""
Refine the bnb box, in particular, 
Given βlb, the maximum tc is bounded by the charging function.
"""
function refine!(box::DsBnBBox, B::Float64)
    b_lb = box.xlb[3] 
    tc_ub = cf_inv(b_lb, B, B)
    box.xub[1] = min(tc_ub, box.xub[1])
    if box.xlb[1] > tc_ub
        box.xlb[1] = tc_ub
    end
    return box
end


"""
divide a box along its longest axis
"""
function divide(box0::DsBnBBox{T}) where T
    # error("")
    # x_diff = box.xub - box.xlb
    # fval, idx = findmax(x_diff)

    # tc, tw, b, tau
    # divide on tw and b have few benefit, we may prefer more idx with 1,4
    # If n_repeated=3, we have [1,2,3,4,1,4,1,4]...
    n_repeated = 1
    div_idx::Int = box0.div_idx
    mod_number = 4+(n_repeated-1)*2
    idx1::Int = mod(div_idx, mod_number) 
    idx::Int = (idx1 == 0) ? mod_number : idx1
    if idx > 4
        idx = isodd(idx) ? 2 : 3
    end
    next_div_idx = box0.div_idx + 1
    if (box0.xlb[idx] == box0.xub[idx])
        idx +=1
        next_div_idx += 1
    end
    box1 = copy(box0)
    box2 = copy(box0)
    mid_pt = (box0.xub[idx] + box0.xlb[idx])/2
    box1.xub[idx] = mid_pt
    box2.xlb[idx] = mid_pt
    box1.div_idx = next_div_idx
    box2.div_idx = next_div_idx
    # xub1::T = T(( (i == idx ? (box0.xub[i]+box0.xlb[i])/2 : box0.xub[i])::Float64 for i in 1:4))
    # xlb2::T = T(( (i == idx ? (box0.xub[i]+box0.xlb[i])/2 : box0.xlb[i])::Float64 for i in 1:4))
    # box1 = DsBnBBox{T}(box0.xlb, xub1, -Inf, Inf, box0.div_idx+1)
    # box2 = DsBnBBox{T}(xlb2, box0.xub, -Inf, Inf, box0.div_idx+1)

    # error("")
    return box1,box2
end

function ds_bnb_update_lbub!(box::DsBnBBox, primal::DsPrimal, dual::DsDual, i::Int, ta, st::DateTime, B, option::DsOption)
    x, box.fub = ds_bnb_get_ub(box, primal, dual, i, ta, st, B, option)
    flb, box.div_idx = ds_bnb_get_lb(box, primal, dual, i, ta, st, B, option)
    ## need to be careful about this!!
    # box.flb = max(flb, box.flb)
    box.flb = flb
    if (box.flb > box.fub + 1e-8)
        @exfiltrate
        error("box.flb > box.fub")
        # @assert (box.flb <= box.fub + 1e-8) "$(box.flb) $(box.fub) $(box.xlb) $(box.xub)"
    end
    return box
end

function ds_bnb_get_ub(box::DsBnBBox, primal::DsPrimal, dual::DsDual, i::Int, ta, st::DateTime, B, option::DsOption)
    #x = (box.xlb + box.xub) /2
    # This is biased, xlb seems likely to achieve lower bound.
    # x = box.xlb
    x_mid = (box.xlb + box.xub) / 2
    val1 = ds_sigma_objective(box.xlb, primal, dual, i, ta, st, B, option)
    val2 = ds_sigma_objective(x_mid, primal, dual, i, ta, st, B, option)
    # val2 = Inf
    x = val1 < val2 ? box.xlb : x_mid
    val = min(val1, val2)
    # val = minimum([ds_sigma_objective(x, primal, dual, i, ta, st, B, option) for x in [box.xlb, box.xub]]) 
    return x, val
    
end

function ds_bnb_get_minmax_rate(tc_lb, tc_ub, b_l, b_u, B)
    b_lb_perc = charge_function(b_l, tc_lb, B) /B
    b_ub_perc = charge_function(b_u, tc_ub, B) /B
    lb_idx = searchsortedlast(g_perc_vec, b_lb_perc)
    ub_idx = searchsortedfirst(g_perc_vec, b_ub_perc)
    min_rate = g_default_rate * g_mul_vec[ub_idx]
    max_rate = g_default_rate * g_mul_vec[lb_idx]
    if (b_ub_perc == 1) 
        # min_rate = 0.0
    end
    if (b_lb_perc == 1) 
        # max_rate = 0.0
    end
    # @show b_lb_perc b_ub_perc lb_idx ub_idx min_rate max_rate
    @assert(min_rate <= max_rate)
    return min_rate, max_rate 
end

function get_price_lbub(tau_l::Float64, tau_u::Float64, ta::TS.TimeArray, st::DateTime)

    π_l::Float64, π_u::Float64 = Inf, -Inf
    # tau = tau_l + tw_l
    # while tau <= tau_u+tw_u
    tau = tau_l
    while tau <= tau_u
        p = get_price(ta, st, tau)[1]
        π_l = min(π_l, p)    
        π_u = max(π_u, p)
        tau += 3600.0
        tau = round(tau/3600.0, RoundDown)*3600.0
    end
    pp = get_price(ta, st, tau_u)[1]
    π_l = min(π_l, pp)    
    π_u = max(π_u, pp)
    
    return π_l, π_u
end

function ds_bnb_get_lb(box::DsBnBBox, primal::DsPrimal, dual::DsDual, i::Int, ta::TS.TimeArray, st::DateTime, B, option::DsOption)
    (;t_mul, b_mul, rho, N) = option
    (tc_l, tw_l, b_l, tau_l) = box.xlb
    (tc_u, tw_u, b_u, tau_u) = box.xub
    (;λ,μ) = dual
    (;beta_cs_vec, tau_cs_vec) = primal
    # price_vec = [get_price(ta, st, tau)[1] for tau in getrange(tau_l+tw_l, tau_u+tw_u, 3600.0)]
    # π_l = minimum(price_vec)
    @assert (tau_l <= tau_u) "$tau_l $tau_u"
    @assert (tw_l <= tw_u) "$tw_l $tw_u"
    π_l::Float64, π_u::Float64 = get_price_lbub(tau_l+tw_l, tau_u+tw_u, ta, st)
    # π_l::Float64, π_u::Float64 = Inf, -Inf
    # tau = tau_l + tw_l
    # while tau <= tau_u+tw_u
    #     p = get_price(ta, st, tau)[1]
    #     π_l = min(π_l, p)    
    #     π_u = max(π_u, p)
    #     tau += 3600.0
    #     tau = round(tau/3600.0, RoundDown)*3600.0
    # end
    # pp = get_price(ta, st, tau_u+tw_u)[1]
    # π_l = min(π_l, pp)    
    # π_u = max(π_u, pp)
    

    # π_l = minimum(tau -> , tau_range)
    # π_u = maximum(tau -> get_price(ta, st, tau)[1], getrange(tau_l+tw_l, tau_u+tw_u, 3600.0))

    λi1 = (i == N+1) ? 0.0 : λ[i+1]
    μi1 = (i == N+1) ? 0.0 : μ[i+1]
    # tc_coeff: we use two linear function to approximate the charging function.
    # println()
    # min_rate, max_rate = ds_bnb_get_minmax_rate(tc_l, tc_u, b_l, b_u, B)
    # println()
    # min_rate_b, max_rate_b = ds_bnb_get_minmax_rate(0.0, 0.0, b_l, b_u, B)
    # ∂ϕ/∂β = (charging rate * 1/(charging rate(no tc))) - 1

    cf_coeff = π_l - μi1
    cf_coeff2 = π_u - μi1
    # if cf_coeff < 0
    #     tc = tc_u; b = b_l
    #     # tc_coeff = cf_coeff * max_rate * b_mul
    # else
    #     # tc_coeff = cf_coeff * min_rate * b_mul
    #     tc = tc_l; b = b_u
    # end

    # tc_coeff = b_mul * min(
    #     cf_coeff * max_rate,
    #     cf_coeff * min_rate,
    # )
    # b_coeff = b_mul * min(
    #     cf_coeff2 * (max_rate/min_rate_b-1),
    #     cf_coeff2 * (min_rate/max_rate_b-1)
    # )
    # a1 = b_mul * cf_coeff * (charge_function(b_l, tc_l, B) - b_l )
    # a2 = - b_coeff * b_l 
    # a3 = - tc_coeff * tc_l

    tc_coeff = b_coeff = 0.0
    a1::Float64 = Inf
    for b in SA[b_l, b_u]
        for tc in SA[tc_l, tc_u]
            val::Float64 = b_mul * cf_coeff * (charge_function(b, tc, B) - b )
            a1 = min(a1, val)
        end
    end
    # vals = (b_mul * cf_coeff * (charge_function(b, tc, B) - b ) for b in (b_l, b_u) for tc in (tc_l, tc_u))
    # a1, a1_idx = findmin(vals)
    # a1 = minimum(vals)
    a2 = 0.0; a3 = 0.0

    # b_mul * cf_coeff * (charge_function(b, tc, B) - b )
    # a2 = t_mul * λi1*( tc_l + tw_l )
    # a3 = t_mul * (λ[i+1]-λ[i]) * ( (λ[i+1] - λ[i] < 0) ? tau_u : tau_l)
    # a4 = b_mul * (μ[i]  -μ[i+1]) * ( (μ[i] - μ[i+1] < 0) ? b_u : b_l)

    last_beta = primal.sub_sol_vec[i].β
    last_tau = primal.sub_sol_vec[i].τ
    beta_cs = max(beta_cs_vec[i+1], 0.0)
    tau_cs = tau_cs_vec[i+1]
    last_tc = primal.sub_sol_vec[i].tc
    last_tw = primal.sub_sol_vec[i].tw

    # augmented panalty to smooth the solution.
    # we want to jointly consider the minimum of the function
    # a1(x-x0)^2 + a2
    a_b = ds_get_augmented_cost(b_l, b_u, beta_cs, rho*b_mul^2, b_mul*(μ[i]-μi1 ) + b_coeff )
    a_tau = ds_get_augmented_cost(tau_l, tau_u, tau_cs, rho*t_mul^2, t_mul*(λi1 - λ[i]))
    a_tc = ds_get_augmented_cost(tc_l, tc_u, last_tc, rho*t_mul^2, t_mul* (λi1) + tc_coeff)
    a_tw = ds_get_augmented_cost(tw_l, tw_u, last_tw, rho*t_mul^2, t_mul* λi1)
    
    # if tau_l <= 

    # c6 = rho * b_mul^2 ( beta_cs_vec[i+1] - β )^2
    # c7 = rho * t_mul^2 ( tau_cs_vec[i+1] - τ )^2

    # @show tc_l, tc_u b_coeff
    # @show min_rate max_rate min_rate_b max_rate_b
    obj_lb = a_b + a_tau + a_tc + a_tw + a1 + a2 + a3
    # _, div_idx = findmin(x->(x), (a_tc, a_tw, a_b, a_tau))
    # x_diff = (box.xub-box.xlb)
    # _, div_idx = findmax(x_diff)
    div_idx = box.div_idx 
    if option.debug
        @show obj_lb a1 a2 a3 a_b a_tau a_tc a_tw cf_coeff b_coeff b_mul*(μ[i]-μi1 ) + b_coeff π_l π_u
        @show beta_cs b_l b_u
        @show (π_u-π_l)/π_u
        # @show dual box.xlb box.xub
    end
    # @show box vals; println()

    return obj_lb::Float64, div_idx::Int
end

"""
minimize the function: a1*(x-x0)^2 + a2*x s.t. x∈[lb,ub]
"""
function ds_get_augmented_cost(lb::T, ub::T, x0::T, a1::T, a2::T) where T <: Float64
    f(x) = a1*(x-x0)^2 + (a2*x)
    flb = f(lb)
    fub = f(ub)
    min_f = min(flb, fub)
    # @show min_f flb fub ""
    (a1 == 0) && return min_f
    # @assert (a1 != 0)
    xs = - a2/2/a1 + x0
    if abs(x0) > 1e6
        # @show xs lb ub x0 a1 a2 flb fub
    end
    if lb <= xs <= ub
        # xs lb ub
        fval = f(xs)
        # @show "test"
        return min(fval, min_f)
    else
        return min_f
    end
    
    # if lb <= x0 <= ub
    #     return 0.0
    # else
    #     c = min( (x0-lb)^2, (x0-ub)^2 )
    #     return multi * c
    # end
        
end