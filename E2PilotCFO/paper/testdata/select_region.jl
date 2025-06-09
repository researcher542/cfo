"""
We want to select a src/des pair in FAF5 for the example such that 
(i) it cross multipl states.
(ii) there are multiple route options.
(iii) the number of nodes is moderate.
(iv) the carbon intensity are largely diverse across the different states.
"""

import E2PilotCore as ep

faf_df = ep.read_faf_data_raw()
faf_reg_dict = ep.read_faf_region()

print_cnt = 0
for (irow, row) in enumerate(eachrow(faf_df))
    src_str = row.dms_orig
    des_str = row.dms_dest
    src_latlon = faf_reg_dict[src_str]
    des_latlon = faf_reg_dict[des_str]
    dist = ep.distance2d(src_latlon, des_latlon)
    if occursin("Rest", src_str) || occursin("Rest", des_str)
        continue
    end
    # if 100e3 <= dist <= 300e3 && occursin("", src_str)
    value = row["million dollars in 2017"]
    # value = row["million ton-miles in 2017"]
    # if occursin("Atlanta", src_str) && value > 1 && dist <= 500e3 && !occursin("GA", des_str)
    # if occursin("AZ", src_str) && value > 1 && dist <= 10000e3 && !occursin("AZ", des_str) && occursin("WA", des_str)
    if value > 1 && dist >= 500 * ep.g_meter_per_mile && dist <= 1000.0 * ep.g_meter_per_mile 
        # && dist >= 500e3 # && !occursin("FL", des_str) && occursin("KY", des_str)
        @info "" irow src_str des_str dist/ep.g_meter_per_mile value
        global print_cnt += 1
    end
    if print_cnt > 10
        break
    end
end

