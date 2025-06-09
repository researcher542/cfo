"""
The data is from the paper 
Synthetic European road freight  transport flow data

Link to paper: https://www.sciencedirect.com/science/article/pii/S235234092101060X?via%3Dihub
Link to dataset: https://data.mendeley.com/datasets/py2zkrb65h/1
"""

using Mmap


function get_etis_data_selected()
    return joinpath(k_data_path, "src.des.pair", "ETIS_OD_selected.csv") 
end

function get_etis_data_zippath()
    return joinpath(k_data_path, "src.des.pair", "ETIS_OD.zip") 
end

function get_etis_raw_df()
    zippath = get_etis_data_zippath()
    archive = ZipReader(mmap(open(zippath)))
    @time "reading csv_file from zip" csv_file = zip_readentry(archive, "01_Trucktrafficflow.csv")
    @time "reading csv file to df" df = CSV.read(csv_file, DataFrame;)
    return df
end

function get_etis_region_df()
    zippath = get_etis_data_zippath()
    archive = ZipReader(mmap(open(zippath)))
    csv_file = zip_readentry(archive, "02_NUTS-3-Regions.csv")
    df = CSV.read(csv_file, DataFrame;)
    return df 
end

function get_etis_region_dict()
    df = get_etis_region_df()
    pair_vec = [] 
    for row in eachrow(df)
        pair = row["ETISPlus_Zone_ID"] => (
            name=row["Name"], 
            lat=row["Geometric_center_Y"], lon=row["Geometric_center_X"], 
            country=row["Country"])
        push!(pair_vec, pair)
    end
    d = Dict(pair_vec)
    return d
end

function getdistance(row, ::ETIS)
    # return row.Total_distance .* 1000.0
    return row.straight_distance 
end

function get_src_latlon(df, ::ETIS)
    lat_vec = df.orig_lat
    lon_vec = df.orig_lon
    return [LatLon(lat, lon) for (lat, lon) in zip(lat_vec, lon_vec)]
end

function get_des_latlon(df, ::ETIS)
    lat_vec = df.dest_lat
    lon_vec = df.dest_lon
    return [LatLon(lat, lon) for (lat, lon) in zip(lat_vec, lon_vec)]
end


function read_od_data_df(::ETIS)
    df = DataFrame(CSV.File(get_etis_data_selected()))
    return df 
end

function get_src_name_vec(df, ::ETIS)
    return df.Name_origin_region
end

function get_des_name_vec(df, ::ETIS)
    return df.Name_destination_region
end

# function read_etis_data_dis_group(net, dis_group::Tuple{T, T},  n_data::Int) where T <: Real
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

