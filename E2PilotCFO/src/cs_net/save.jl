


function load_cs_net(region::String, cs_nei_flag::Bool; tmp=false, ow_flag::Bool = false)
    path = get_cs_net_path(region, cs_nei_flag, tmp) 
    @info "loading cs_net $path"
    default_cs_net = CsNet(-1; cs_nei_flag = cs_nei_flag)
    if !isfile(path)
        @warn "cs_net not found: $path. Loading the empty cs_net. You may need to run the get_cs_net.jl to get the cs_net."
    end
    if !isfile(path) || ow_flag
        return default_cs_net
    end
    #data_dict = load(g_bin_cs_net_path)
    data_dict = JLD2.load(path)
    if !haskey(data_dict, region)
        return default_cs_net
    end
    return data_dict[region]
end

"""
Get the file paht of stored cs_net
"""
function get_cs_net_path(region::String, cs_nei_flag::Bool, tmp::Bool)
    bin_cs_net_path = joinpath(k_data_path, "map", "bin", "cs_net")
    path = bin_cs_net_path * "-" * region
    if cs_nei_flag
        path = path * "-nei"
    end
    if tmp
        path = path * "-tmp"
    end
    path *= ".jld2"
    return path
end

let
timer = now()
cnt = 0

global function save_cs_net(region::String, cs_net::CsNet, cs_nei_flag::Bool; tmp=true)
    path = get_cs_net_path(region, cs_nei_flag, tmp) 
    @info "saving cs_net tmp=$(tmp), path=$path please do not exit."
    if isfile(path)
        data_dict = JLD2.load(path) 
        data_dict[region] = cs_net
    else
        data_dict = Dict(region => cs_net)
    end
    if (now() - timer >= Minute(30)) && cnt >= 10
        @debug "backing up cs_net"
        cp(path, path*"$(now()).backup"; force=true)
        timer = now()
        cnt = 0
    end
    cnt += 1
    JLD2.save(path, data_dict)
    @info "cs_net saved."
end

end