"""
Functions to analyze and manipulate the data from the origin-destination Pairs
"""

abstract type AbstractODSet end

struct FAF <: AbstractODSet end
struct ETIS <: AbstractODSet end



include("fafdata.jl")
include("etisdata.jl")

"""
Read OD data from a specific distance group
"""
function read_od_data_dis_group(net, dis_group::Tuple{T, T},  n_data::Int, odset::AbstractODSet) where T <: Real
    # longhaul_distance = ep.g_longhaul_distance
    min_dis, max_dis = dis_group
    df = read_od_data_df(odset)
    df1 = filter(
        row -> getdistance(row, odset) > min_dis * g_meter_per_mile && 
        getdistance(row, odset) <= max_dis * g_meter_per_mile, 
        df)

    if odset == FAF()
        df1 = sort!(df1, "million dollars in 2017", rev=true)
    elseif odset == ETIS()
        df1 = sort!(df1, "Traffic_flow_tons_2019", rev=true)
    end

    src_vec = get_src_latlon(df1, odset)
    des_vec = get_des_latlon(df1, odset)
    src_idx_vec = closest_node(net, src_vec)
    des_idx_vec = closest_node(net, des_vec)
    src_name_vec = get_src_name_vec(df1, odset)
    des_name_vec = get_des_name_vec(df1, odset)
    src_des_vec = [
        (src=src_idx_vec[i], 
         des=des_idx_vec[i],
         src_name = src_name_vec[i],
         des_name = des_name_vec[i],
         ) for i in 1:n_data
    ]

    return src_des_vec
end