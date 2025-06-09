"""
Utils to handle OD datasets
"""

g_dis_group_vec = [
    (500.0, 1000.0),
    (1000.0, 1500.0),
    (1500.0, 2000.0),
    (2000.0, 8000.0)
]

g_dis_group_emission_vec = [
    # 65.32 # 250-500
    40.49+25.05, # 500-750 and 750-1000
    31.9, # 1000-1500
    16.22, # 1500-2000
    25.74, # > 2000
]

g_dis_group_emission_dict = OrderedDict(
    g_dis_group_vec .=> g_dis_group_emission_vec
)

function disgroup2label(min_dis, max_dis, odset = FAF())
    dis_group_vec = get_dis_group_vec(odset)
    if min_dis == dis_group_vec[end][1]
        # label = @sprintf "\$\\geq  %d\$" min_dis
        label = @sprintf "%d+" min_dis
    else
        label = @sprintf "%d-%d" min_dis max_dis
    end
    return label
    
end

function get_dis_group_vec(::FAF)
    return g_dis_group_vec
end

function get_dis_group_vec(::ETIS)
    return [
        (500.0, 1000.0),
        (1000.0, 1500.0),
        (1500.0, 8000.0),
    ]
end


function read_od_data_all_group(net, n_data_per_group::Int, odset::AbstractODSet) 
    src_des_vec = []
    for dis_group in get_dis_group_vec(odset)
        src_des_vec1 = ep.read_od_data_dis_group(net, dis_group, n_data_per_group, odset)
        src_des_vec = append!(src_des_vec, src_des_vec1)
    end
    return collect(src_des_vec)
end

# function read_faf_data_all_group(net, n_data_per_group::Int) 
#     src_des_vec = Tuple{Int, Int}[]
#     for dis_group in g_dis_group_vec
#         src_des_vec1 = read_faf_data_dis_group(net, dis_group, n_data_per_group)
#         src_des_vec = append!(src_des_vec, src_des_vec1)
#     end
#     return src_des_vec
# end

# function read_faf_data_dis_group(net, dis_group::Tuple{T, T},  n_data::Int) where T <: Real
#     # longhaul_distance = ep.g_longhaul_distance
#     min_dis, max_dis = dis_group
#     df = ep.read_faf_data_df()
#     df1 = filter(
#         row -> row.distance > min_dis * ep.g_meter_per_mile && 
#         row.distance <= max_dis * ep.g_meter_per_mile, 
#         df)
# 
#     df1 = sort!(df1, "million dollars in 2017", rev=true)
# 
#     src_vec = df1.orig_latlon
#     des_vec = df1.dest_latlon
#     src_idx_vec = closest_node(net, src_vec)
#     des_idx_vec = closest_node(net, des_vec)
#     src_des_vec = [
#         (src=src_idx_vec[i], 
#          des=des_idx_vec[i],
#          src_name = df1.dms_orig[i],
#          des_name = df1.dms_dest[i],
#          ) for i in 1:n_data
#     ]
# 
#     return src_des_vec
# end

