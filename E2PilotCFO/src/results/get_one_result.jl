
"""
T_flag: if true, use the specified deadline.
"""
function get_one_result!(prob0::AbsCfoProb, alg_name::String, res0::CfoResult; t_ratio::Real=1.2, fast_res = nothing,  option0, ice_veh = nothing)
    # @info "get_one_result!" alg_name
    prob = copy(prob0)
    res = deepcopy(res0)
    option = deepcopy(option0)
    (;src, des) = res
    if src == des
        @warn "Get the same src and des. Simply return." src des
        return res
    end
    prob.src = src
    prob.des = des
    res.alg_name = alg_name
    if hasfield(typeof(option), :beta_lb)
        res.beta_lb = option.β_lb
    end
    res.region = prob.net.region
    if occursin("fast", alg_name)
        prob.objtype = ObjTime()
    elseif occursin("ds-e", alg_name)
        prob.objtype = ObjEnergy()
    else
        prob.objtype = ObjCarbon()
    end
    # prob.objtype = occursin("ds-e", alg_name) ? ep.ObjEnergy() : ep.ObjCarbon()
    t_begin = time()
    if alg_name == "fast"
        primal = ds_fast_path(prob, option)
    else
        if !isnothing(fast_res)
            prob.T = fast_res.time * t_ratio 
            @info "Using fast_res.time=$(fast_res.time) to set the deadline T=$(prob.T)"
            if occursin("ds", alg_name)
                option.primal0 = fast_res.primal
            end
        end
        T = prob.T
        res.deadline = T
        res.t_ratio = t_ratio
        if occursin("ds", alg_name)
            state_vec = cfo_dual_subgradient(prob, option)
            # kp_option = cKpOption(;K=8, restrict_N = false, fast_flag=true, thread_flag=true)
            # state_vec = cfo_kpath(prob, option, kp_option)
            # state_vec = cfo_bundle_method(prob, option)
            state = state_vec[end]
            primal = state_vec[end].primal_f
            if !isvalid(primal)
                @warn "ds failed to get a valid solution, revert to fastest result."
                primal = deepcopy(fast_res.primal)
            end
            res.lb = state.theta_c
        elseif occursin("practice", alg_name)
            primal, _ = cfo_practice(prob)
        elseif occursin("reopt", alg_name)
            ## If we want to reoptimize the the fast_res
            @assert !isnothing(fast_res)
            if occursin("reopt-spd-wait", alg_name) 
                @debug alg_name "Computing reopt-spd-wait"
                dp_option = DsReoptDpOption(;speed_opt_flag=true, charge_opt_flag=true)
                _, primal = ds_reopt_dp(prob, fast_res.primal, option, dp_option)
                # _, primal = ds_reopt(prob, fast_res.primal, option)
            elseif occursin("reopt-spd", alg_name) 
                @debug alg_name "Computing reopt-spd"
                dp_option = DsReoptDpOption(;speed_opt_flag=true, charge_opt_flag=false)
                _, primal = ds_reopt_dp(prob, fast_res.primal, option, dp_option)
                # _, primal = ds_reopt(prob, fast_res.primal, option)
            end
        elseif occursin("tbexnet", alg_name)
            stats1 = @timed exp_net = construct_bte_net(prob, option)
            res.btexnet_construct_time = stats1.time
            res.option = option
            res.exp_nv = length(exp_net.nodesdata)
            res.exp_ne = length(exp_net.w_dict)
            res.orig_nv = length(prob.net.nodesdata)
            res.orig_ne = length(prob.net.waydata)
            try
                stats2 = @timed primal = cfo_time_battery_expanded(prob, option, exp_net)
                res.btexnet_solve_time = stats2.time
            catch err
                err_msg = sprint(showerror, err)
                res.error = err_msg
                return res
            end
        elseif occursin("grb", alg_name)
            res.orig_nv = length(prob.net.nodesdata)
            res.orig_ne = length(prob.net.waydata)
            primal, model = cfo_mip_jump(prob, option)
            res.model_status = termination_status(model)
        elseif occursin("ice", alg_name)
        # if alg_name == "paso-ice"
            (;net,src,des) = prob
            if alg_name == "fast-ice"
                paso_step = ep.paso(net, src, des, 0; veh=ice_veh, breakearly=true)
            else
                paso_step = ep.paso(net, src, des, T; veh=ice_veh, breakearly=true)
            end
            res.time = paso_step.summary.duration
            res.cost = paso_step.summary.cost
            res.step = paso_step
            return res
        else
            error("unknown algorithm name: $(alg_name)!")
        end
    end
    # path_new, t_vec_new, β_vec, τ_vec = ep.polish(net, path, t_vec, ev, β0, start_time; γ=γ)
    #obj1, β_vec1, τ_vec1 = ep.simulate(net, path_new, t_vec_new, ev, β0, start_time, 1e20; check_flag=true)
    prob.objtype = ep.ObjCarbon()

    sim_res = ds_simulate(prob, primal; predict_mode=PredictPerfect())
    (;obj, beta_vec, tau_vec, infeasible_flag, e_cost) = sim_res
    ################
    # obj, β_vec1, τ_vec1, t_vec1 = ep.simulate_cs_net(prob, path, t_vec)
    res.time = tau_vec[end] / option.t_mul
    res.cost = obj / option.b_mul
    res.tau_vec = tau_vec
    res.β_vec = beta_vec / option.b_mul
    res.primal = primal
    res.infeasible_flag = infeasible_flag
    res.e_cost = e_cost
    res.cpu_time = time() - t_begin ## in seconds
    res.carbon_dataset = prob.carbon_dataset
    res.predict_mode = prob.predict_mode
    res.start_time = prob.start_time
    # res.t_vec = t_vec1
    # res.path = path
    # res.ori_time = sum(t_vec)
    return res 
end

cfo_mip_jump() = nothing

"""
Find the OD idx of a given od_vec

region has the format of "src-des"
"""
function find_od_idx(od_vec, region::String)
    parts = split(region, "-")
    src_idx = parse(Int, parts[1])
    des_idx = parse(Int, parts[2])
    # @show parts
    src_str = @sprintf "%03d" src_idx
    des_str = @sprintf "%03d" des_idx
    @show src_str, des_str
    function is_target_od(od)
        flag1 = occursin(src_str, od.src_name)
        flag2 = occursin(des_str, od.des_name)
        return flag1 && flag2
    end
    idx = findfirst(is_target_od, od_vec)
    return idx
end

"""
Get the result for one instance, dedicated for the method comparison.

In this 
"""
function get_one_result_cmp(prob0::AbsCfoProb, prob_us::AbsCfoProb, option, alg_name0::String; t_ratio::Real = 1.2, ow_flag::Bool = false, fast_res = nothing)

    prob = copy(prob0)
    if occursin("ds", alg_name0)
        # prob = copy(prob_us)
    end
    (; src, des, start_time, β0, ev) = prob

    t_begin = time()
    # alg_name0 = "tbexnet"
    res0 = CfoResult(; alg=alg_name0, src=src, des=des, idx=0, start_time=start_time, cap=ev.cap, β0=β0)
    res0.region = prob0.net.region
    res0.option = option
    res0.od_data = FAF()
    res0.orig_nv = length(prob.net.nodesdata)
    res0.orig_ne = length(prob.net.waydata)

    file_name = get_result_file_name(res0)
    file_path = joinpath(g_cfo_result_dir, file_name)

    local res_new = deepcopy(res0)
    if (!ow_flag) && isfile(file_path)
        @info "skip $file_path"
        return JLD2.load_object(file_path)
    end

    @info "solving $file_path with ow_flag=$ow_flag, isfile=$(isfile(file_path))"
    try 
        res_new = get_one_result!(prob, alg_name0, res0; t_ratio = t_ratio, option0=option, fast_res=fast_res) 
    catch err
        # res_new.error = "$err"
        if isa(err, InterruptException)
            @warn "Got interrupt exception" err
            rethrow(err)
            return res_new
        end
        err_msg = sprint(showerror, err)
        res_new.error = err_msg
        res_new.cpu_time = time() - t_begin
        bt = catch_backtrace()
    	# res_new.backtrace = bt
        res_new.region = prob0.net.region
        @error "Got error with" prob0 res0.idx prob0.start_time res0
    	showerror(stdout, err, bt)
        # rethrow(err)
    end
    @info "saving result to $file_path"
    res_new.region = prob0.net.region
    save_one_result(res_new) 
    flush(stdout)
    GC.gc(true)
    return res_new
    
end