"""
Combine the sperate results to the compact result.
Some command line functions to manupilate the results files.
"""

function get_combined_file_name(region::String)
    combined_file_name = "cfo-results.$(region).$(today()).jld2"
    return combined_file_name
end

"""
cp_flag: if true, copy the combined result to the result path.
"""
function combine_cfo_results(region::String, save::Bool = false, cp_flag::Bool = false)
    @info "Combining cfo results for $region."
    res_dict = OrderedDict{String, Vector{CfoResult}}()
    combined_file_name = get_combined_file_name(region)
    nfile = length(readdir(g_cfo_result_dir; join=true))
    @show nfile combined_file_name
    all_res_vec = fill(CfoResult(), nfile)
    @sync for (ifile, file_path) in enumerate(readdir(g_cfo_result_dir; join=true))
        if ifile % 100 == 0
            @info "reading result $(ifile)/$(nfile)"
        end
        if occursin(combined_file_name, file_path)
            # If it is the combined file.
            continue
        end
        if !occursin(region, file_path)
            continue
        end
        if occursin("cfo-results", file_path)
            # If it is the cfo result.
            @warn "Skip the cfo result file $file_path."
            continue
        end
        Threads.@spawn begin
            # @info "Loading result file $file_path."
            res = JLD2.load_object(file_path)
            if isnothing(res)
                @warn "Error in loading $file_path."
            end
            all_res_vec[ifile] = res
        end
        # res = try 
        #     res = JLD2.load_object(file_path)
        # catch e
        #     @warn "Error in loading $file_path." e
        #     continue
        # end
        
    end
    @show length(filter(x->!isempty(x.alg), all_res_vec))

    for res in all_res_vec
        res_key = res.alg
        if isempty(res.alg)
            continue
        end
        if haskey(res.meta, :renewable_mul) && res.renewable_mul != 1.0
            # ratio = renewable_multiplier_to_ratio(res.renewable_mul)
            # res_key *= @sprintf "-renew-ratio=%.3f" res.renewable_mul
            res_key *= @sprintf "-mul=%.3f" res.renewable_mul
        end
        if !haskey(res_dict, res_key)
            res_dict[res_key] = CfoResult[]
        end
        push!(res_dict[res_key], res)
    end

    for (key, res_vec) in res_dict
        res_vec1 = sort!(res_vec, by=res->res.idx)
        res_dict[key] = res_vec1
    end

    if save
        save_path = joinpath(g_cfo_result_dir, combined_file_name)
        JLD2.save_object(save_path, res_dict)
        if cp_flag
            result_dir = joinpath(ep.k_data_path, "results")
            dst_path =  joinpath(result_dir, combined_file_name)
            @info "copying the combined result to $dst_path."
            cp(save_path, dst_path, force=true)
        end
        @info "NOTE: We also need to change the g_cfo_result_path to $dst_path."
    end
    return res_dict
end

function seprate_cfo_results(;tmp=true, verbose::Bool=true, overwrite::Bool=true)
    res_dict = load_cfo_result(;tmp=tmp)
    n = sum(length(res_vec) for res_vec in values(res_dict))
    @info "Seperating $n results."
    idx = 0
    for (key, res_vec) in res_dict
        for res in res_vec
            idx += 1
            if verbose
                @info "Seperating $idx/$n result."
            end
            save_one_result(res, overwrite)
        end
    end
end

