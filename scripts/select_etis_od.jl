"""
A script to select ETIS OD data based on a tonnage threshold and save the selected data to a CSV file.
"""

import E2PilotCore as ep
import E2PilotCFO as cfo
# import E2PilotServer as es
using CSV, HTTP, DataFrames, Geodesy, Statistics, Graphs, PyPlot, ZipArchives, Mmap
# e2plt = Base.get_extension(cfo, :PlotExt)
# e2plt.set_pyplot_style_no_latex()
pygui(false)


if !isdefined(Main, :df)
    df = ep.get_etis_raw_df()
end
if !isdefined(Main, :prob)
    prob = cfo.get_test_prob("eu")
end

@time region_dict = ep.get_etis_region_dict()
### The box of EU regions by OSM, we can run the scripts/osm2csv.jl to get those bounding box

(;min_lat, max_lat, min_lon, max_lon) = ep.get_eu_osm_box()


region_id_set = Set(
    reg_id for (reg_id, reg) in region_dict if (
        reg.lat >= min_lat && reg.lat <= max_lat && 
        reg.lon >= min_lon && reg.lon <= max_lon
    )
)



df1 = filter(row -> row["Total_distance"] > 100.0, df)
df2 = filter(row -> row["ID_origin_region"] in region_id_set && 
    row["ID_destination_region"] in region_id_set, df1)
col_idx = 2
col = ["Traffic_flow_tons_2019", "Traffic_flow_trucks_2019"][col_idx]; 
val_tol = [500.0, 50.0][col_idx]


df3 = filter(row -> row[col] > val_tol, df2)
@time df4 = sort(df3, col, rev=true)
deleted_cols = ["Edge_path_E_road", "Unnamed: 0"]
for col in deleted_cols
    select!(df4, Not(col))
end

selected_path = ep.get_etis_data_selected()


# @time region_df = ep.get_etis_region_df()

df4.orig_lat = [region_dict[row."ID_origin_region"].lat for row in eachrow(df4)]
df4.orig_lon = [region_dict[row."ID_origin_region"].lon for row in eachrow(df4)]
df4.dest_lat = [region_dict[row."ID_destination_region"].lat for row in eachrow(df4)]
df4.dest_lon = [region_dict[row."ID_destination_region"].lon for row in eachrow(df4)]

function get_distance_net(net, row)
    latlon_o = ep.LatLon(row.orig_lat, row.orig_lon)
    nd1 = ep.closest_node(net, latlon_o)

    latlon_d = ep.LatLon(row.dest_lat, row.dest_lon)
    nd2 = ep.closest_node(net, latlon_d) 

    dis = ep.distance2d(net, nd1, nd2) / ep.g_meter_per_mile
    return dis
end

function get_distance_net_vec(net, df)
    latlon_o_vec = [ep.LatLon(row.orig_lat, row.orig_lon) for row in eachrow(df)]
    latlon_d_vec = [ep.LatLon(row.dest_lat, row.dest_lon) for row in eachrow(df)]

    o_vec = ep.closest_node(net, latlon_o_vec)
    d_vec = ep.closest_node(net, latlon_d_vec)

    dis_vec = ep.distance2d.((net,), o_vec, d_vec)
    return dis_vec
end

@time df4.straight_distance = get_distance_net_vec(prob.net, df4)

@time df4 = filter(row -> row.straight_distance > 500.0 * ep.g_meter_per_mile, df4)

# for row in eachrow(df4)
    # @time dis = get_distance_net(prob.net, row)
    # @show dis 
# end

# df4.straight_distance = [ep.distance2d(LatLon(row.orig_lat, row.orig_lon), LatLon(row.dest_lat, row.dest_lon)) for row in eachrow(df4)]
# @time df4.straight_distance = [get_distance_net(prob, row) for row in eachrow(df4)]

# for row in eachrow(df4)
#     dis1 = ep.distance2d(
#         LatLon(row.orig_lat, row.orig_lon), 
#         LatLon(row.dest_lat, row.dest_lon))
#     dis2 = row.Total_distance * 1000.0
# end

df5 = filter(row -> ep.getdistance(row, ep.ETIS()) > 1500 * ep.g_meter_per_mile, df4)
@show length(df4.Total_distance)
@show length(df5.Total_distance)

CSV.write(selected_path, df4)

@info "done."


