"""
remove part of results that are produced by some buggy codes.
"""


function remove_all_results(region::String, dry_run::Bool = true)
    @info "removing all results with dry_run=$dry_run."
    function check_func(res)
        return true
    end 
    return _remove_results(region, check_func, dry_run)
end

"""
"""
function remove_practice_results(region::String, dry_run::Bool=true)
    @info "removing all results which are generated from practice method."
    function check_func(res)
        return occursin("practice", res.alg)
    end

    return _remove_results(region, check_func, dry_run)
end

"""
"""
function remove_reopt_results(region::String, dry_run::Bool=true)
    @info "removing all results which are generated from reopt method."
    function check_func(res)
        return occursin("reopt", res.alg)
    end

    return _remove_results(region, check_func, dry_run)
end

function remove_all_ds_results(region::String, dry_run::Bool = true)
    @info "removing all results which are generated from dual subgradient method."
    function check_func(res)
        flag = occursin("ds", res.alg)
        # @show flag res.alg
        return flag
    end 
    return _remove_results(region, check_func, dry_run)
end

function remove_infeasible_results(region::String, dry_run::Bool = true)
    @info "removing results which are infeasible."
    function check_func(res)
        if occursin("ds", res.alg)
            return res.infeasible_flag
        else
            return false
        end
    end 
    return _remove_results(region, check_func, dry_run)
end

function is_alpha_consistent(res::CfoResult)
    alg_str = res.alg
    if !occursin("alpha", alg_str)
        return true
    end
    alpha_str = alg_str[end-4:end]
    alpha = parse(Float64, alpha_str)
    if !isvalid(res)
        @warn "Got a invalid res" res.alg res.idx
        @exfiltrate
        error("")
    end
    alpha0 = res.beta_lb / res.cap
    if !(alpha ≈ alpha0)
        return false
        # @debug res.alg alpha alpha0 res.idx
    end
    return true
end

"""
Due to some bugs, some old result files do not have a modified β_lb. We need to re-run those instances.
"""
function remove_inconsistent_alpha_results(region::String, dry_run::Bool=true)
    @info "removing results with in-consistent results with β_lb."
    function check_func(res)
        return !is_alpha_consistent(res) 
    end

    return _remove_results(region, check_func, dry_run)
end



"""
region: test4node, TXOH, oldmap
"""
function remove_invalid_results(region::String, dry_run::Bool=true)
    @info "removing results with error."
    function check_func(res)
        return !isvalid(res) 
    end
    return _remove_results(region, check_func, dry_run)
end

"""
check_func: CfoResult -> Bool
If check_func is true, remove the cooresponding result file.
This is a low_level API
"""
function _remove_results(region::String, check_func, dry_run::Bool = true)
    combined_file_name = get_combined_file_name(region)
    rm_file_vec = String[]
    for (ifile, file_path) in enumerate(readdir(g_cfo_result_dir; join=true))
        if occursin(combined_file_name, file_path) || (!occursin("jld2", file_path)) 
            # If it is the combined file.
            continue
        end
        if (!occursin("idx", file_path))
            continue
        end
        if !occursin(region, file_path)
            continue
        end
        # @show ifile file_path
        res = JLD2.load_object(file_path)
        if check_func(res)
            push!(rm_file_vec, file_path)
        end
    end
    @warn "removing $(length(rm_file_vec)) files. dry_run=$(dry_run)"
    
    if !dry_run
        for filepath in rm_file_vec
            rm(filepath)
        end
    else
        for f in rm_file_vec
            @show f
        end
        return rm_file_vec
    end
    
end