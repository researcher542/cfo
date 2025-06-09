using Distributed


"""
Get the results for multiple instances.
"""
function get_results(prob::AbsCfoProb, n_data::Int, alg_name::String ; ow_flag=false, t_ratio=1.2, option, ice_veh = nothing, alpha::Float64 = 0.05, nthread::Int= 1)
    GC.gc()
    src_des_vec = cfo.read_faf_data_all_group(net, n_data)
    res_dict = load_cfo_result(; tmp=true)
    option.β_lb = alpha * prob.ev.cap
    
    alg_name = get_result_key(alg_name, t_ratio, alpha)
    res_vec = get_exist_result(alg_name, ow_flag, res_dict) 
    max_task = nthread*1
    arg_vec = []

    function _get_one_result(args)
        (prob0, alg_name0, res0, t_ratio0, option0, fast_res0, ice_veh0) = args
        local res_new = deepcopy(res0)
        try 
            res_new = get_one_result!(prob0, alg_name0, res0; t_ratio = t_ratio0, option0=option0, fast_res=fast_res0, ice_veh = ice_veh0) 
        catch err
            res_new.error = err
        end
        
        return res_new
    end
    
    for i = 1:n_data
        src,des = src_des_vec[i]
        src,des = closest_cs_node.((net,), (src,des))
        res = CfoResult(; alg=alg_name, src=src, des=des, idx=i, start_time=start_time, cap=ev.cap, β0=prob.β0)
        if (!ow_flag) && already_has(res, res_vec)
            @info "skip $alg_name for i=$i alpha=$alpha"
            continue
        end
        @info "computing for $alg_name path i=$i" src des 
        if alg_name == "fast" || alg_name == "fast-ice"
            fast_res = nothing
        elseif occursin("fast-ds", alg_name)
            fast_res = res_dict["fast"][i]
        else
            fast_res = res_dict["fast-ds"][i]
        end
        
        
        # res_task = Distributed.@spawnat :any get_one_result!(prob, alg_name, res; t_ratio = t_ratio, option=option, fast_res=fast_res, ice_veh = ice_veh) 

        push!(arg_vec, (prob, alg_name, deepcopy(res), t_ratio, deepcopy(option), deepcopy(fast_res), ice_veh))
        if length(arg_vec) >= max_task
            @debug "Fetching and waiting for results." length(arg_vec)
            #res_vec1 = fetch.(arg_vec)
            # res_vec1 = ep.e2map(_get_one_result, arg_vec)
            res_vec1 = pmap(_get_one_result, arg_vec)
            res_vec = vcat(res_vec, res_vec1)
            # idx_vec = [res.idx for res in res_vec]
            # @assert all(idx_vec .== collect(1:length(idx_vec)))
            arg_vec = []
            res_dict[alg_name] = res_vec
            save_cfo_result(res_dict; tmp=true)
            @info "Results got with $(length(res_vec))/$(n_data)."
        end
    end
    res_vec1 = pmap(_get_one_result, arg_vec)
    res_vec = vcat(res_vec, res_vec1)
    res_dict[alg_name] = res_vec
    if !isempty(arg_vec)
        save_cfo_result(res_dict; tmp=true)
    end
end

function already_has(res0::CfoResult, res_vec::Vector)
    if haskey(res0.meta, :error)
        return false
    end
    for res in res_vec 
        if res == res0
            return true
        end
    end
    return false
end

"""
Get the results for multiple instances.

data_idx_vec: the index of the data we want to get for each start time.
"""
function get_exist_result_all_st(alg_name::String, t_ratio::Float64, data_idx_vec; kwargs...)
    all_vec = []
    st_year = g_default_st_year
    if haskey(kwargs, :st_year)
        st_year = kwargs[:st_year]
    end

    for st in g_st_vec
        st1 = ep.change_year(st, st_year)
        res_vec = get_exist_result(alg_name, t_ratio, st1; kwargs...)
        append!(all_vec, res_vec[data_idx_vec]) 
        # @show st, alg_name, length(res_vec[data_idx_vec])
        # @show length(all_vec)
    end
    return all_vec
end

"""
The high-level function to get the results for multiple instances.
"""
function get_exist_result(alg_name::String, t_ratio::Float64, st::DateTime,
    ;
    ow_flag::Bool = false, 
    res_dict = nothing, 
    alpha::Float64 = 0.05, 
    renewable_mul::Float64 = 1.0,
    scenario_key = g_default_scenario_key, 
    st_year = g_default_st_year, 
    set_infeasible = true, 
    predict_mode = PredictPerfect(), 
    carbon_dataset = CambiumDataset(),
    prob = nothing,
    )
    st0 = deepcopy(st)
    ## The name key in the res_dict...
    if !is_ds_carbon_result_key(alg_name) && (predict_mode != PredictNoise())
        st = g_default_st
    end
    if st_year != g_default_st_year
        st = ep.change_year(st, st_year)
    end
    # if !is_ice_result_key(alg_name)
    #     # st = DateTime(st.year, st.month, st.day, 0, 0, 0)
    # end

    name_key = get_result_key(alg_name, t_ratio, alpha, st; 
        scenario=scenario_key, st_year=st_year, predict_mode=predict_mode, carbon_dataset=carbon_dataset
        )
    
    res_vec0 = get_exist_result(name_key, ow_flag, res_dict, renewable_mul; set_infeasible=set_infeasible)
    res_vec = deepcopy(res_vec0)

    ## If it is not carbon optimized operation for e-truck, we need to recompute the result by its start time. This is because we skip the cases of energy-optimized solution in computing the results to save computational time.
    # if !is_ds_carbon_result_key(alg_name) && (st0 != g_default_st) && 
    if is_etruck_result_key(alg_name)
    # if true
        # @info "Recompute the result for $alg_name with st=$st0"
        # @show alg_name
        if isnothing(prob)
            @error("prob is nothing, please provide a prob.")
        end
        
        for (ires, res) in enumerate(res_vec)
            if ismissing(res)
                continue
            end
            prob1 = copy(prob)
            prob1.src = res.src
            prob1.des = res.des
            if !occursin("fast", alg_name)
                prob1.T = res.deadline
            end
            prob1.start_time = st0
            # @show scenario_key 
            prob1.scenario = g_scenario_dict[scenario_key]
            prob1.predict_mode = predict_mode
            # prob1.carbon_dataset = carbon_dataset

            sim_res = ds_simulate(prob1, res.primal; check_flag = false, predict_mode = prob1.predict_mode)
            # @show res.start_time, st0
            # @show res.cost, sim_res.obj * 3.6e6
            res_vec[ires].cost = sim_res.obj * 3.6e6
            res_vec[ires].start_time = st0
            # @show sim_res
            # (ires > 3) && break
        end
    end
    return res_vec
end


function is_ds_carbon_result_key(k::AbstractString)
    return occursin("ds-c", k) 
end

function is_ice_result_key(k::AbstractString)
    return occursin("-ice", k) && !occursin("practice", k)
end

function is_etruck_result_key(k::AbstractString)
    return !is_ice_result_key(k)
end

"""
Get the results for multiple instances.

set_infeaible: if result is infeasible, we set the cost to be NaN and do not compare it.
"""
function get_exist_result(
    name0::String, ow_flag::Bool, res_dict = nothing, 
    renewable_mul::Float64 = 1.0
    ; 
    set_infeasible = true)
    name = deepcopy(name0)
    if renewable_mul != 1.0
        name *= @sprintf "-mul=%.3f" renewable_mul
        @show name
    end
    if ow_flag 
        res_vec = Vector{CfoResult}()
    else
        if isnothing(res_dict)
            res_dict = load_cfo_result(;tmp=false)
        end
        if haskey(res_dict, name)
            res_vec = res_dict[name]
        else
            @warn "res_dict does not contain key name=$name."
            res_vec = Vector{CfoResult}()
            return res_vec
        end
    end
    res_vec0 = sort(res_vec, by=res -> res.idx)
    n = res_vec0[end].idx
    res_vec = Vector{Any}(missing, n)
    for res in res_vec0
        if res.infeasible_flag && set_infeasible
            # @show res.idx
            res1 = deepcopy(res)
            res1.cost = NaN
            res1.e_cost = NaN
            res_vec[res.idx] = res1
        else
            res_vec[res.idx] = res
        end
    end
    return res_vec
end

function get_result_key(alg_name::String, t_ratio::Float64, alpha::Float64, st::DateTime
    ; 
    scenario = g_default_scenario_key, 
    st_year = g_default_st_year, 
    predict_mode = PredictPerfect(), 
    carbon_dataset = CambiumDataset())
    if !occursin("fast", alg_name) && !occursin("practice", alg_name)
        alg_name = @sprintf "%s-ratio=%.3f-st=%s" alg_name t_ratio st
        if alpha != 0.05
            alg_name = @sprintf "%s-alpha=%.3f" alg_name alpha
        end
        if scenario != g_default_scenario_key
            alg_name = @sprintf "%s-scenario=%s" alg_name scenario
        end
        if st_year != g_default_st_year
            alg_name = @sprintf "%s-st_year=%d" alg_name st_year
        end
        if carbon_dataset != CambiumDataset()
            alg_name = @sprintf "%s-%s" alg_name carbon_dataset2str(carbon_dataset)
            if predict_mode != PredictPerfect()
                alg_name = @sprintf "%s-predict=%s" alg_name prediction_mode2str(predict_mode)
            end
        end
    end
    return alg_name
end

"""
Previously I forgot to add the e_cost, need to append it here.
"""
function update_result_energy_cost!(prob::AbsCfoProb, res_dict)
    src0 = prob.src
    des0 = prob.des
    β0 = prob.β0
    for (key,res_vec) in res_dict
        if occursin("ice", key) 
            continue
        end
        for res in res_vec
            if haskey(res, :e_cost)
                continue
            end
            prob.src = res.src
            prob.des = res.des
            prob.β0 = res.β0
            (;e_cost) = ep.ds_simulate(prob, res.primal; check_flag = false, predict_mode = prob.predict_mode)
            res.e_cost = e_cost
        end
    end
    prob.src = src0
    prob.des = des0
    prob.β0 = prob.β0
    return res_dict
end



