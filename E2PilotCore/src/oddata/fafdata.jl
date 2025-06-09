"""
utility to process the Freight Analysis Framework (FAF) data
We use the following website to get the faf data
# https://faf.ornl.gov/faf5/dtt_domestic.aspx
"""

const g_faf_data_dir = joinpath(k_data_path, "src.des.pair")

"""
read the raw data from the FAF website.
"""
function read_faf_data_raw()
    faf_data_path = joinpath(g_faf_data_dir, "faf_selected_subset.csv")
	df = DataFrame(CSV.File(faf_data_path))
    df = sort!(df, "million dollars in 2017", rev=true)
    return df 
end



"""
"""
function read_faf_data_df()
    df = read_faf_data_raw()
    reg_dict = read_faf_region()
    orig_latlon_vec = [reg_dict[dms_orig] for dms_orig in df.dms_orig]
    df[!, "orig_latlon"] = orig_latlon_vec

    dest_latlon_vec = [reg_dict[dms_dest] for dms_dest in df.dms_dest]
    df[!, "dest_latlon"] = dest_latlon_vec

    distance_vec = [distance2d(row.orig_latlon, row.dest_latlon) for row in eachrow(df)]
    df[!, "distance"] = distance_vec

    return df
end

function read_faf_data_longhaul(n_sample::Int = 1000, longhaul_distance::Float64 = g_longhaul_distance)
    df = read_faf_data_df()
    df1 = filter(row -> row.distance > longhaul_distance, df)
    df2 = df1[1:n_sample, :]
    faf_data_vec = collect(zip(df2.orig_latlon, df2.dest_latlon, df2."million dollars in 2017"))
    return faf_data_vec
end

"""
read the data converted by the python script.
"""
function read_faf_data_old(n_sample::Int = 1000)
    faf_data_path = joinpath(g_faf_data_dir, "src_des_pair.csv")
	df = DataFrame(CSV.File(faf_data_path))
    faf_data_vec = []
    for (irow, row) in enumerate(eachrow(df))
        if irow > n_sample break end
        latlon_src = LatLon(row.src_lat, row.src_lon)
        latlon_des = LatLon(row.des_lat, row.des_lon)
        push!(faf_data_vec, [latlon_src, latlon_des, row.value])
    end
    return faf_data_vec
end

function read_faf_data_net(net, n_data::Int)
    # faf_data_vec = read_faf_data_longhaul(n_data)
    faf_data_vec = read_faf_data_old(n_data)
    src_vec = [data[1] for data in faf_data_vec]
    des_vec = [data[2] for data in faf_data_vec]
    src_idx_vec = closest_node(net, src_vec)
    des_idx_vec = closest_node(net, des_vec)
    src_des_vec = [
        (src_idx_vec[i], des_idx_vec[i]) for i in 1:n_data
    ]
    # dis_vec = [
    #     distance2d(net, src_idx_vec[i], des_idx_vec[i]) for i in 1:n_data
    # ]
    # sort!(dis_vec)
    return src_des_vec
end

let

local_faf_region_dict::Dict{String, LatLon{Float64}} = Dict{String, LatLon{Float64}}()
"""
read the FAF5 region data, return a dict map from region to latlon
"""
global function read_faf_region()

    if isempty(local_faf_region_dict)
        faf_data_path = joinpath(g_faf_data_dir, "faf_reg_latlon.csv")
	    df = DataFrame(CSV.File(faf_data_path))
        for (irow, row) in enumerate(eachrow(df))
           latlon = LatLon(row.lat, row.lon)
           local_faf_region_dict[row.region] = latlon
           # latlon_src = LatLon(row.src_lat, row.src_lon)
           # latlon_des = LatLon(row.des_lat, row.des_lon)
           # push!(faf_data_vec, [latlon_src, latlon_des, row.value])
       end
    end
    return local_faf_region_dict
end
end

function read_df_region_vec()
    dict = read_faf_region()
    return keys(dict)
end

function fafregion2idx(s::AbstractString)
    parts = split(s, "-")  # Split the string at the "-"
    number = parse(Int, parts[1])  # Convert the first part to an integer
    return number
end

function idx2fafregion(idx0::Int)
    for reg in read_df_region_vec()
        idx = fafregion2idx(reg)
        if idx0 == idx
            return reg
        end
    end
    return ""
end

function fafidx2netidx(net, idx::Int)
    faf_region = idx2fafregion(idx) 
    region_dict = read_faf_region()
    latlon = region_dict[faf_region]
    net_idx = closest_node(net, latlon)
    return net_idx
end

function get_src_latlon(df, ::FAF)
    return df.orig_latlon 
end

function get_des_latlon(df, ::FAF)
    return df.dest_latlon 
end

function get_src_name_vec(df, ::FAF)
    return df.dms_orig
end

function get_des_name_vec(df, ::FAF)
    return df.dms_dest
end

function read_od_data_df(::FAF)
    return read_faf_data_df() 
end

function getdistance(row, ::FAF)
    return row.distance
end