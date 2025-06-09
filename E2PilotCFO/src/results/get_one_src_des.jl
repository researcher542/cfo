


function update_faf_src_des!(prob::AbsCfoProb, faf_idx::Int)
    net = prob.net
    src_des_vec = read_faf_data_all_group(net, 100)
    res = src_des_vec[faf_idx] 
    src0 = res.src
    des0 = res.des
    src, des = closest_cs_node.((net,), (src0, des0))
    prob.src = src
    prob.des = des
    return src, des
end

function read_od_data_vec(prob::AbsCfoProb)
    net = prob.net
    od_vec = read_od_data_all_group(net, 100, prob.odset)
    # if occursin("eu", region)
    #     # src_des_vec = read_faf_data_dis_group(net, (500.0, 1000.0), 100)
    # else
    #     src_des_vec = read_faf_data_all_group(net, 100) 
    # end
    return od_vec
end

function read_od_data(prob::AbsCfoProb, od_idx::Int)
    src_des_vec = read_od_data_vec(prob)
    res = src_des_vec[od_idx] 
    return res
end

"""
Get all the results of one source destination pair and save it to the result directory with idx.
"""
function get_one_src_des(prob0::AbsCfoProb, idx::Int, ice_veh, option0; ow_flag::Bool = false, dry_run::Bool = false, run_alpha::Bool = false, run_t_ratio::Bool = false, run_future::Bool = false, run_noise::Bool = false)
    prob = copy(prob0)
    @info "Getting results for idx=$idx region=$(prob.net.region)" 
    GC.gc(false)
    option = deepcopy(option0)
    (;net, ev) = prob
    start_time = g_default_st
    prob.start_time = start_time
    if net.region == "test4node"
        od_data = nothing
        src, des = 1, 4
    else
        od_data = read_od_data(prob, idx)
        src = closest_cs_node(net, od_data.src)
        des = closest_cs_node(net, od_data.des)
        
        # src, des = update_faf_src_des!(prob, idx)
        # src_des_vec = read_faf_data_all_group(net, 100)
        # src0, des0 = src_des_vec[idx]
        # src, des = closest_cs_node.((net,), (src0, des0))
    end

    
    res = CfoResult(; alg="", src=src, des=des, idx=idx, start_time=start_time, cap=ev.cap, β0=prob.β0)
    res.region = prob.net.region
    res.od_data = od_data
    arg_vec = []
    g_alpha = 0.05
    g_t_ratio = 1.2
    for alg_name in ["fast",  "fast-ice", "practice"]
        fast_res = nothing
        args = (prob, alg_name, deepcopy(res), g_t_ratio, option, fast_res, ice_veh, idx, ow_flag)
        push!(arg_vec, args)
    end
    res_vec = ep.e2map(_get_one_result, arg_vec)

    alg_name = "fast-ds"
    fast_res = res_vec[1]
    args = (prob, alg_name, deepcopy(res), g_t_ratio, option, fast_res, ice_veh, idx, ow_flag)
    fast_ds_res = _get_one_result(args)

    

    arg_vec = []
    ### The base case
    for alg in ["ds-c", "ds-e", "paso-ice"]
        prob1 = copy(prob)
        st_vec = [g_default_st]
        if alg == "ds-c"
            st_vec = g_st_vec
        end
        for st in st_vec
            prob1.start_time = st
            alg_name = get_result_key(alg, g_t_ratio, g_alpha, prob1.start_time)
            args = (prob1, alg_name, deepcopy(res), g_t_ratio, deepcopy(option), fast_ds_res, ice_veh, idx, ow_flag)
            push!(arg_vec, args)
        end
    end

    if run_alpha
        for alpha in g_alpha_vec
            option = deepcopy(option0)
            option.β_lb = alpha * prob.ev.cap
            alg = "ds-c"
            alg_name = get_result_key(alg, g_t_ratio, alpha, start_time)
            args = (prob, alg_name, deepcopy(res), g_t_ratio, deepcopy(option), fast_ds_res, ice_veh, idx, ow_flag)
            push!(arg_vec, args)
        end
    end

    if run_t_ratio
        option.β_lb = g_alpha * prob.ev.cap
        for t_ratio in g_t_ratio_vec
            # for alg in ["ds-c", "ds-e", "paso-ice", "reopt-spd", "reopt-spd-wait"]
            for alg in ["ds-c", "ds-e", "paso-ice"]
                if alg == "ds-c"
                    for st in g_st_vec
                        prob1 = copy(prob)
                        prob1.start_time = st
                        alg_name = get_result_key(alg, t_ratio, g_alpha, prob1.start_time)
                        args = (prob1, alg_name, deepcopy(res), t_ratio, deepcopy(option), fast_ds_res, ice_veh, idx, ow_flag)
                        push!(arg_vec, args)
                    end
                else
                    prob1 = copy(prob)
                    prob1.start_time = g_default_st
                    alg_name = get_result_key(alg, t_ratio, g_alpha, prob1.start_time)
                    args = (prob1, alg_name, deepcopy(res), t_ratio, deepcopy(option), fast_ds_res, ice_veh, idx, ow_flag)
                    push!(arg_vec, args)
                end
            end
        end
    end

    if run_future
        option.β_lb = g_alpha * prob.ev.cap
        # for (key, scenario) in g_scenario_dict
        for scenario_key in g_selected_scenario_vec
            scenario = g_scenario_dict[scenario_key]
            for st_year in g_st_year_vec
                for st in g_st_vec
                    for alg in ["ds-c"]
                        option = deepcopy(option0)
                        prob1 = copy(prob)
                        prob1.scenario = scenario
                        # st = prob1.start_time

                        prob1.start_time = DateTime(st_year,
                            Month(st).value, Day(st).value, 
                            Hour(st).value, Minute(st).value, Second(st).value)

                        alg_name = get_result_key(alg, g_t_ratio, g_alpha, prob1.start_time; scenario=scenario_key, st_year=st_year)
                        args = (prob1, alg_name, deepcopy(res), g_t_ratio, deepcopy(option), fast_ds_res, ice_veh, idx, ow_flag)
                        push!(arg_vec, args)
                    end
                end
            end
        end
    end

    

    # We need to be careful with this part. since the fix data also needs to be changed.
    if run_noise
        carbon_dataset = CarbonCastDataset()
        prob1 = copy(prob)
        prob1.start_time = DateTime(2022, 8, 1, 12, 0, 0)
        res1 = deepcopy(res)
        res1.start_time = prob1.start_time
        net1 = deepcopy(prob1.net)
        empty!(net1.meta)
        net1.sep_cs_node_flag = net.sep_cs_node_flag
        # cs_vec = get_charge_station(net1, prob1.scenario; fuel=false, carbon_dataset=carbon_dataset)
        cs_vec1 = deepcopy(net.cs_vec)
        change_carbon_dataset!(cs_vec1, carbon_dataset, prob1.scenario, )
       
        add_charge_station!(net1, cs_vec1)
        cs_nei_flag = false
        cs_net = load_cs_net(net1.region, cs_nei_flag)
        cs_net = cs_net_restrict_cap!(cs_net, ev.cap)
        net1.cs_net = cs_net
        net1.continent = prob.net.continent

        
        noise_fix_data = CfoProbData(
           net=net1, ev=ev, 
           carbon_dataset=carbon_dataset, 
           carbon_dict= get_carbon_dict(carbon_dataset),
           carbon_predict_dict = get_carbon_dict_predict(carbon_dataset),
           cs_dict = net1.cs_dict
        )

        prob1.fix_data = Ref(noise_fix_data)
        t_ratio = g_t_ratio
        for pred_mode in [PredictPerfect(), PredictNoise()]
            alg_name = get_result_key("ds-c", t_ratio, g_alpha, prob1.start_time; predict_mode=pred_mode, carbon_dataset = carbon_dataset)
            prob1.predict_mode = pred_mode
            args = (deepcopy(prob1), alg_name, deepcopy(res1), t_ratio, deepcopy(option), fast_ds_res, ice_veh, idx, ow_flag)
            push!(arg_vec, args)
        end
    end

    @info "Got $(length(arg_vec)) instances"
    # nthread = round(Int, Threads.nthreads() / 32)
    # nthread = max(nthread, 1)
    # nthread = 1
    # @info "starting e2map with nthread=$nthread..."

    if dry_run
        return arg_vec
    else
        for (iarg, arg) in enumerate(arg_vec)
            @info "In get_one_src_des, computing iarg=$iarg of $(length(arg_vec))"
            _get_one_result(arg)
            GC.gc(true)
        end
        # res_vec = ep.e2map(_get_one_result, arg_vec, true, true; prefix = "get_one_src_des (idx=$idx): ", nthread=nthread)
    end
     
    @info "Got results for idx=$idx." 

    return res_vec
end

"""
"""
function _get_one_result(args)
    # @debug "solving args" args[2]
    (prob0, alg_name0, res0, t_ratio0, option0, fast_res0, ice_veh0, idx, ow_flag) = args
    # file_name = get_result_file_name(idx, alg_name0, prob0.net.region)
    res0.alg = alg_name0

    file_name = get_result_file_name(res0)
    file_path = joinpath(g_cfo_result_dir, file_name)
    local res_new = deepcopy(res0)
    if (!ow_flag) && isfile(file_path)
        @info "skip $file_path"
        return JLD2.load_object(file_path)
    else
        @info "solving $file_path"
    end
    try 
        res_new = get_one_result!(prob0, alg_name0, res0; t_ratio = t_ratio0, option0=option0, fast_res=fast_res0, ice_veh = ice_veh0) 
    catch err
        res_new.error = err
        bt = catch_backtrace()
    	# res_new.backtrace = bt
        res_new.region = prob0.net.region
        if !isnothing(fast_res0)
            res_new.deadline = fast_res0.time * t_ratio0 
        end
        @error "Got error with" prob0 res0.idx prob0.start_time res0
    	showerror(stdout, err, bt)
        rethrow(err)
    end
    if !is_alpha_consistent(res_new)
        msg = "res has inconsistent alpha."
        @warn msg
        error(msg)
    end
    # @debug "res_new" args[2] res_new
    @info "saving idx=$idx $alg_name0."
    save_one_result(res_new) 
    @info "idx=$idx $alg_name0 saved."
    flush(stdout)
    GC.gc(true)
    return res_new
end