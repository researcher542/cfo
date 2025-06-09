"""
The utility functions to save and load results etc.
"""

include("config.jl")
include("resultdata.jl")
include("get_one_result.jl")
include("get_results.jl")
include("convert_results.jl")
include("loadsave.jl")
include("get_one_src_des.jl")
include("remove_results.jl")


"""
If the current environment good for plotting.
"""
function is_plot_env()
    if ep.is_in_hpc()
        return false
    else
        return true
    end
end

function get_src_des(net, faf_data_vec, i)
    (latlon_src,latlon_des) = faf_data_vec[i]
    src = ep.closest_node_in_junction(net, latlon_src)
    des = ep.closest_node_in_junction(net, latlon_des)
    return (src, des)
end

function get_fast_cfo_result(idx::Int, src, des::Int)
    res_vec1 = load(g_cfo_result_path, "fast")
    res = res_vec1[idx]
    @assert src == res.src 
    @assert des == res.des
    return res
end
