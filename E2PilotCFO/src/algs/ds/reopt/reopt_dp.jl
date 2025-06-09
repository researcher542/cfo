"""
Use DP to reoptimize the problem.
"""

@with_kw struct DsReoptDpState
    F::Float64
    beta::Float64
    tau::Float64
    istage::Int
    inflag::Bool # If true, it represent the state before charging.
end

# @with_kw struct DsReoptScoreMat
# end

function Base.lt(o::Base.ForwardOrdering, s1::DsReoptDpState, s2::DsReoptDpState)
    return (s1.F < s2.F) ||
           ((s1.F == s2.F) && (s1.istage < s2.istage)) ||
           ((s1.F == s2.F) && (s1.istage == s2.istage) && (s1.tau < s2.tau)) ||
           ((s1.F == s2.F) && (s1.istage == s2.istage) && (s1.tau == s2.tau) && (s1.beta > s2.beta))
end

function get_ibeta_itau(s::DsReoptDpState, eps_b::Float64, eps_t::Float64)
    ibeta::Int = round(Int, s.beta / eps_b, RoundDown)
    itau::Int = round(Int, s.tau / eps_t, RoundUp)
    return ibeta,itau 
end

function Base.:(==)(s1::T, s2::T) where T <: DsReoptDpState
    return  (s1.F == s2.F) && (s1.beta == s2.beta) && (s1.tau == s2.tau) && (s1.istage == s2.istage) && (s1.inflag == s2.inflag)
    # return  (s1.beta == s2.beta) && (s1.tau == s2.tau) && (s1.istage == s2.istage) && (s1.inflag == s2.inflag)
end

"""
propogate the state for the road segment
"""
function ds_get_next_state_tr(s::DsReoptDpState, tr::Float64, e_cost::Float64, dF::Float64)
    @assert (!s.inflag)
    return DsReoptDpState(
        s.F + dF, s.beta - e_cost, s.tau + tr, s.istage, true
    )
end

"""
propogate the state for charging station
    db: Δβ
"""
function ds_get_next_state_cs(s::DsReoptDpState, tw::Float64, tc::Float64, db::Float64, dF::Float64)
    @assert (s.inflag)
    return DsReoptDpState(
        s.F + dF, s.beta + db, s.tau + tw + tc, s.istage + 1, false
    )
end

function getkey(s::DsReoptDpState)
    return (s.beta, s.tau, s.istage, s.inflag) 
end

"""
Push the state for the charging edge.
"""
function push_state!(pq, score_mat, pre_state_dict, cur_state, next_state, eps_b, eps_t, B, T, b_lb, debug)

    beta_flag = (0 <= next_state.beta <= B)
    tau_flag = (0 <= next_state.tau <= T)
    blb_flag = (next_state.beta < b_lb) && (cur_state.inflag == false)

    if debug && cur_state.istage == 2 && cur_state.inflag == false
        # @warn "pushing state" cur_state next_state beta_flag tau_flag blb_flag 
    end
        
    if (!beta_flag) || (!tau_flag) || (blb_flag)
        return
    end
    ibeta_, itau_ = get_ibeta_itau(next_state, eps_b, eps_t)
    istage_ = next_state.istage
    # score_key = (istage_, ibeta_, itau_, next_state.inflag)
    # score_ = get(score_mat, score_key, Inf)
    inflag_idx = (next_state.inflag ? 1 : 0)
    score_key =  CartesianIndex(istage_, ibeta_, itau_,  inflag_idx)
    score_ = score_mat[score_key]
    score_flag = (next_state.F < score_)

       
    if score_flag
        score_mat[score_key] = next_state.F
        pre_state_dict[next_state] = cur_state
        # enqueue!(pq, next_state, next_state)
        # pq[next_state] = next_state
        pq[getkey(next_state)] = next_state
    end
end

"""
push the states for the traveling edge.
"""
function push_state_travel!(prob, nstage::Int, cs_net, u, v, objtype, pq, score_mat, pre_state_dict, cur_state, eps_b, eps_t, B, T, b_lb, debug, speed_opt_flag, init_ci::Real)

    # cs_edge::CsEdge{Float64} = get_edge(cs_net, cs_path[istage], cs_path[istage+1])
    cs_edge::CsEdge{Float64} = get_edge(cs_net, u, v)
    (;min_t, max_t) = cs_edge
    i_min_t = round(Int, min_t / eps_t, RoundUp)
    i_max_t = round(Int, max_t / eps_t, RoundDown) 
    i_max_t = max(i_max_t, i_min_t) # i_min_t > i_max_t can happen when eps_t is not small enough.
    # debug && @debug "" istage min_t max_t i_min_t i_max_t cs_edge
    if (!speed_opt_flag)
        # If we do not optimize speed, simply choose the fastest one.
        i_max_t = i_min_t
    end
    for itr in i_min_t:i_max_t
        tr = itr * eps_t
        e_cost0::Float64 = cs_edge(tr, true)
        e_cost::Float64 = round(Int, e_cost0 / eps_b, RoundUp) * eps_b
        dF = 0.0
        if objtype == ObjTime()
            dF = tr
        else # if objtype == ObjCarbon()
            ## If this is the next stage is the final stage, we need to consider the battery cost at the beginning.
            if cur_state.istage == nstage
                dF = get_dF_final_state(prob, cur_state, e_cost, init_ci) 
                # @show init_delta_beta ci dF
            end
        end
        next_state = ds_get_next_state_tr(cur_state, tr, e_cost, dF)

        push_state!(pq, score_mat, pre_state_dict, cur_state, next_state, eps_b, eps_t, B, T, b_lb, debug)
        # push_state!(cur_state, next_state)
    end # end of for loop for travel time
    
end

function push_state_charge!(prob::AbsCfoProb, Nt, Nb, cs_idx::Int, cur_itau, cur_ibeta, pq, score_mat, pre_state_dict, cur_state::DsReoptDpState, eps_b, eps_t, B, T, b_lb, debug, charge_opt_flag)
    for ib in cur_ibeta:Nb
        next_b = ib * eps_b
        tc = cf_inv(cur_state.beta, next_b, B)

        if (!charge_opt_flag)
            # If we do not optimize charge, simply wait for minimum amount of time possible.
            tau = cur_state.tau + tc + g_min_wait_time
            itau0 = round(Int, tau/eps_t, RoundUp)
            itau_range = itau0:itau0
        else
            itau_range = cur_itau:Nt
        end
        # @debug "" charge_opt_flag cur_itau Nt itau_range
        cfo_obj_state = CfoObjState()
        for itau::Int in itau_range
            next_tau::Float64 = itau * eps_t
            tw::Float64 = next_tau - cur_state.tau - tc
            if (tw < g_min_wait_time) || (tw > g_max_wait_time)
                continue
            end
            db = next_b - cur_state.beta
            # dF = cfo_objective(prob, cs_idx, cur_state.beta, tc, cur_state.tau + tw, prob.predict_mode, cfo_obj_state)
            predict_mode = prob.predict_mode
            dF = cfo_objective(prob, cs_idx, cur_state.beta, tc, cur_state.tau + tw, predict_mode, cfo_obj_state)
            if prob.objtype == ObjTime()
                ## Note that the objective is the time, the cfo_objective already handles the charging time
                dF += tw
            end
            # @show itau prob.objtype tw
            if tw > g_min_wait_time
                next_state = ds_get_next_state_cs(cur_state, tw, tc, db, dF) 
                push_state!(pq, score_mat, pre_state_dict, cur_state, next_state, eps_b, eps_t, B, T, b_lb, debug)
                # push_state!(cur_state, next_state)
            end
        end
    end
    # @assert false
end

@with_kw struct DsReoptDpOption
    Nb::Int = 100
    eps_t::Float64 = 300.0
    speed_opt_flag::Bool = true # If we want to do the speed optimization, if false, we always choose the fastest path.
    charge_opt_flag::Bool = true # If we want to do the charing optimization, if false, we do not wait.
    timeout::Float64 = 10 * 60.0 # in seconds. By default, we set it to 10 minutes.
end

"""
Re-optimize primal solution with self-implemented dynamic programming.
Nb: the number of states of battery.
eps_t: the time step in seconds.
"""
function ds_reopt_dp(prob::AbsCfoProb, primal::DsPrimal, option::DsOption, dp_option::DsReoptDpOption = DsReoptDpOption(); debug::Bool = false)
    (;Nb, eps_t, speed_opt_flag, charge_opt_flag, timeout) = dp_option
    (;src, des, β0, T, net) = prob
    (;cs_net) = net
    (;b_mul) = option

    t_begin = time()

    B = prob.ev.cap
    Nt = round(Int, T/eps_t, RoundDown)
    eps_b = B / Nb
    beta0 = round(Int, β0 / eps_b, RoundUp)* eps_b
    cs_path = get_cs_path(primal)
    nstage = length(cs_path) - 1
    # pq = HeapPriorityQueue{DsReoptDpState, DsReoptDpState}()
    pq = PriorityQueue{Tuple{Float64, Float64, Int, Bool}, DsReoptDpState}()
    state0 = DsReoptDpState(0.0, beta0, 0.0, 1, false)
    # enqueue!(pq, state0, state0)
    pq[getkey(state0)] = state0
    des_state = DsReoptDpState(0.0, 0.0, 0.0, -1, false)
    ib_lb = round(Int, option.β_lb/eps_b, RoundUp)
    b_lb = ib_lb*eps_b
    if option.verbose > 0
        # @debug "In ds_reopt_dp" Nt Nb nstage
    end
    max_size = (Nt*Nb*nstage*2)

    # store the minimum score of current stage
    ## istage, ibeta, itau, inflag
    # score_mat = OrderedDict{Tuple{Int, Int, Int, Bool}, Float64}() # fill(Inf, nstage, Nb, Nt)
    score_mat = OA.Origin(0)(fill(Inf, nstage+1, Nb+1, Nt+2, 2))
    # the dict to store the previous information.
    pre_state_dict = OrderedDict{DsReoptDpState, DsReoptDpState}()
    # sizehint!(score_mat, max_size)
    sizehint!(pre_state_dict, max_size)

    cnt::Int = 0
    max_n_state::Int = 0
    max_mem_usage::Float64 = (Sys.total_memory() - Sys.free_memory())/1e6
    init_ci = get_init_battery_ci(prob, prob.predict_mode)
    while !isempty(pq)
        # cur_state = dequeue!(pq)
        _, cur_state = dequeue_pair!(pq)
        cnt += 1
        if (cur_state.istage == nstage)  && (cur_state.inflag == true)
            des_state = cur_state
            # @show cur_state
            # @debug "break"
            break
        end
        if (time() - t_begin) > timeout
            @warn "Timeout in reopt_dp" timeout
            break
        end

        if (cnt % 10000 == 0) && debug
            @debug "" cnt length(pq) Nb Nt nstage
            @show cur_state
            max_n_state = max(max_n_state, length(pq))
            max_mem_usage = max(max_mem_usage, (Sys.total_memory() - Sys.free_memory())/1e6)
            @show max_n_state max_mem_usage
        end
        if (cnt % 100_000 == 0) && debug
            # break 
        end
        
        (;istage)  = cur_state
        cur_itau = round(Int, cur_state.tau/eps_t) 
        cur_ibeta = round(Int, cur_state.beta/eps_b) 
        if cur_state.inflag
            # If at charging station
            if cs_path[istage] == cs_path[istage+1]
                # @debug "same cs."
                tw = tc = db = dF = 0.0
                next_state = ds_get_next_state_cs(cur_state, tw, tc, db, dF) 
                # push_state!(cur_state, next_state)
                push_state!(pq, score_mat, pre_state_dict, cur_state, next_state, eps_b, eps_t, B, T, b_lb, debug)
                continue
            end
            # @debug "" cur_itau cur_ibeta
            # propogate the next state

            cs_idx = cs_path[istage+1]
            push_state_charge!(prob, Nt, Nb, cs_idx, cur_itau, cur_ibeta, pq, score_mat, pre_state_dict, cur_state, eps_b, eps_t, B, T, b_lb, debug, charge_opt_flag)
        else
            # If needs to decide the travel time
            if cs_path[istage] == cs_path[istage+1]
                ## if the next stage is the same charging station, we skip it and go to the next stage.
                tr = dF = e_cost = 0.0
                if cur_state.istage == nstage 
                    dF = get_dF_final_state(prob, cur_state, e_cost, init_ci)
                end
                next_state = ds_get_next_state_tr(cur_state, tr, e_cost, dF) 
                push_state!(pq, score_mat, pre_state_dict, cur_state, next_state, eps_b, eps_t, B, T, b_lb, debug)
                # push_state!(cur_state, next_state)
                continue
            end

            u = cs_path[istage]
            v = cs_path[istage+1]
            push_state_travel!(prob, nstage, cs_net, u, v, prob.objtype, pq, score_mat, pre_state_dict, cur_state, eps_b, eps_t, B, T, b_lb, debug, speed_opt_flag, init_ci)
            
        end # end of if statement for inflag
    end # end of while loop

    if option.verbose > 0
        # @debug "In reopt_dp, max $(max_n_state) states in pq. max_mem_usage $(max_mem_usage) MB" Nt Nb nstage (Nt*Nb*nstage*2) 
        # @show Base.summarysize(pre_state_dict)/1e6 Base.summarysize(score_mat)/1e6
    end

    if des_state.istage == -1
        if option.verbose > 0
            @warn "Cannot find a feasible path in reopt_dp. cs_path=$cs_path" maxlog=10
        end
        return Inf, primal
    end
    ### start to recover the solution.
    cost = des_state.F 
    if prob.objtype != ObjTime()
        cost = cost * option.b_mul
    else
        cost = cost * option.t_mul
    end
    state_vec = DsReoptDpState[des_state]
    cur_state = des_state
    while true
        cur_state = pre_state_dict[cur_state]
        pushfirst!(state_vec, cur_state)
        if (cur_state.istage == 1) && (cur_state.inflag == false)
            break
        end
    end

    t_vec = Float64[]
    tw_vec = Float64[0.0]
    tc_vec = Float64[0.0]

    for istate in 1:length(state_vec)-1
        s1 = state_vec[istate]
        s2 = state_vec[istate+1]
        @assert xor(s1.inflag, s2.inflag)
        dF = (s2.F - s1.F) * b_mul
        if s1.inflag
            tc = cf_inv(s1.beta, s2.beta, B)
            tw = s2.tau - s1.tau - tc
            push!(tw_vec, tw)
            push!(tc_vec, tc)
        else
            tr = s2.tau - s1.tau
            push!(t_vec, tr)
        end
    end

    t_all_vec = vcat(t_vec, tc_vec, tw_vec)
    # @show t_vec tc_vec tw_vec

    ## Mannually enforce GC
    # empty!(score_mat)
    empty!(pre_state_dict)
    empty!(pq)
    # GC.safepoint()
    GC.gc()
    
    # @timeit g_to "tvec2primal" 
    stat = @timed primal_new = ds_reopt_cs_net_tvec2primal(t_all_vec, primal, prob, option)
    if option.verbose > 0
        # print_timed_stat(stat)
    end
    if cost == 0.0
        @warn "cost is zero in reopt_dp" primal_new
    end
    return cost, primal_new
end

function get_dF_final_state(prob::AbsCfoProb, cur_state, e_cost, init_ci)
    final_beta = cur_state.beta - e_cost
    init_delta_beta = prob.β0 - final_beta
    dF = 0.0
    if prob.objtype == ObjTime()
        return 0.0
    elseif prob.objtype == ObjCarbon()
        # ta = get_price_ta(prob, prob.src, prob.predict_mode)
        # init_ci = get_price(ta, prob.start_time, 0.0)[1]
        # init_ci = get_init_battery_ci(prob, prob.predict_mode)
        dF = init_ci * init_delta_beta
    elseif prob.objtype == ObjEnergy()
        dF = init_delta_beta
    end
    
    return dF
end