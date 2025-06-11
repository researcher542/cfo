


"""
Get the number of completed instances.
"""
function get_num_of_done()
    idx::Int = 0
    while true
        idx += 1
        file_name = @sprintf "idx%d.reopt-spd-wait-ratio=1.500.jld2" idx
        file_path = joinpath(g_cfo_result_dir, file_name)
        if !isfile(file_path)
            break
        end
    end
    return idx
end

prediction_mode2str(predict_mode::PredictPerfect) = "perfect"
prediction_mode2str(predict_mode::PredictNoise) = "noise"
prediction_mode2str(predict_mode::PredictML) = "learning"
carbon_dataset2str(carbon_dataset::CambiumDataset) = "cambium"
carbon_dataset2str(carbon_dataset::CarbonCastDataset) = "cast"


# function _get_result_file_name(idx::Int, alg_name::String, region::String, renewable_multiplier::Float64, predict_mode::AbstractPredictionMode)
#     return file_name
# end

function get_result_file_name(res::CfoResult)
    (;idx, alg, region, renewable_mul) = res
    alg_name = alg
    file_name = @sprintf "idx%d.%s.%s" idx alg_name region
    
    if renewable_mul != 1.0
        file_name *= @sprintf ".mul=%.3f" renewable_mul
    end

    file_name *= ".jld2"
    return file_name
end

"""
"""
function save_one_result(res::CfoResult, overwrite::Bool = true; folder::String = "")
    file_name = get_result_file_name(res)
    if folder == ""
        folder = g_cfo_result_dir
    end
        
    file_path = joinpath(folder, file_name)
    if isfile(file_path) && !overwrite
        @info "The file $file_path already exists. Skip it."
        return
    end
    JLD2.save_object(file_path, res)
end

function load_cfo_result(path = nothing; tmp::Bool = true)
    if isnothing(path)
        path = g_cfo_result_path
    end
    if tmp
        path =  path * ".tmp.jld2"
    end
    if isfile(path)
        res_dict = JLD2.load(path)
        if haskey(res_dict, "single_stored_object")
            return res_dict["single_stored_object"]
        end
    else
        res_dict = Dict()
    end
    return res_dict
end

function save_cfo_result(res_dict, path = nothing; tmp::Bool = true)
    if isnothing(path)
        path = g_cfo_result_path
    end
    if tmp
        path =  path * ".tmp.jld2"
    end
    @debug "Saving cfo results: $path."
    JLD2.save(path, res_dict)
    @debug "cfo result saved."
end

function remove_result_file()
    @warn "removing the result file for cfo task. press enter to continue..."
    readline()
    rm(g_cfo_result_path, force=true) 
    rm(g_cfo_result_path * ".tmp.jld2", force=true) 
end