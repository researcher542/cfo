
# ENV["JULIA_DEBUG"] = "Main"

# t_ratio_vec = Vector(1.2:0.1:1.4)

n_data = 100
ow_flag = false
ow_fast_flag = false
n_worker = 5
multi_flag = true
remove_flag = false # remove the whole result file and run the simulation again.
# g_cfo_result_path = joinpath(ep.k_data_path, "results", "cfo-results.jld2") 
global g_result_continent::Symbol = [:us, :eu][1]
g_cfo_result_dir = joinpath(ep.k_data_path, "results", "cfo") 
g_cfo_cmp_result_dir = joinpath(ep.k_data_path, "results", "cfo_cmp") 
# global g_cfo_result_path = joinpath(ep.k_data_path, "results", "cfo-results.dis500.jld2") 
global g_cfo_result_path = joinpath(ep.k_data_path, "results", "cfo-results.oldmap.2025-04-21.jld2") 
fig_dir = ep.k_fig_dir
g_alg_names = ["fast", "paso-ice", "fast-ice", "ds-e", "ds-c", "reopt-spd", "reopt-spd-wait", "fast-ds", "practice"]
g_alg_xlabels = ["PRACTICE_OLD", "PASO-ICE", "FAST-ICE", "ENERGY", "CARBON", "FAST-SPD", "FAST-SPD-WAIT", "FAST", "PRACTICE"]
g_t_ratio_vec = Vector(1.1:0.1:1.5)
g_alpha_vec = Vector(0:0.02:0.12)
g_st_year_vec = [2024, 2026, 2028, 2030, 2035, 2040, 2045, 2050]

## test for different start time
g_st_vec = [
    DateTime(2024, 2, 1, 8, 0, 0),
    DateTime(2024, 5, 1, 8, 0, 0),
    DateTime(2024, 8, 1, 8, 0, 0),
    DateTime(2024, 11, 1, 8, 0, 0),
]
g_default_st = g_st_vec[1]

# g_renewable_ratio_vec = Vector([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
g_renewable_ratio_vec = Vector([0.2, 0.4, 0.6, 0.8])

################ Test for network with 4 nodes.
function get_cmp_od_region_vec(prob,)
    @assert prob.net.region == "oldmap"

    n_data_per_group::Int = 20

    # od_data_vec = read_od_data_all_group(prob.net, n_data_per_group, prob.odset)
    odset = prob.odset
    od_data_vec = []
    # for dis_group in get_dis_group_vec(odset)
    # for dis_group in [(500.0, 1000.0), (1000.0, 1500.0), (1500.0, 2000.0)]
    # for dis_group in [(500.0, 1000.0), (1000.0, 1500.0)]
    # for dis_group in [(500.0, 1000.0), (1000.0, 1500.0)]
    #     od_data_vec = append!(od_data_vec, src_des_vec1)
    # end

    if true
        ## the ODs used in the paper.
        od_region_vec = [
            "171-179" 
            "061-069"  
            "486-484"  
            "121-122"  
            "061-041"  

            "559-551"  
            "064-411"  
            "484-485"  
            "486-489"  
            "124-129"  

            "132-131"  
            "489-486"  
            "171-559"  
            "129-122"  
            # "061-081"  

            "261-452"  
            "081-171"  
            "131-122"  
            "081-089"  
            "292-421"  

            "292-241" # 6539 nodes
            # "486-350" # 5920 nodes
        ]
        return od_region_vec
    else
        od_vec1 = ep.read_od_data_dis_group(prob.net, (100.0, 300.0), 10, odset)
        od_vec2 = ep.read_od_data_dis_group(prob.net, (300.0, 500.0), 10, odset)
        od_vec3 = ep.read_od_data_dis_group(prob.net, (500.0, 1000.0), 10, odset)
        od_vec4 = ep.read_od_data_dis_group(prob.net, (800.0, 1000.0), 10, odset)
        od_vec5 = ep.read_od_data_dis_group(prob.net, (600.0, 800.0), 10, odset)
        # od_vec4 = ep.read_od_data_dis_group(prob.net, (1000.0, 1500.0), 10, odset)
        od_data_vec = vcat(od_vec1, od_vec2, od_vec3, od_vec4, od_vec5)
    end

    unique!(od_data_vec) # remove duplicates

    function oddata2str(od)
        s1 = od.src_name[1:3] 
        s2 = od.des_name[1:3]
        return "$s1-$s2"
    end

    od_region_vec = [
        oddata2str(od) for od in od_data_vec
    ]
    # od_region_vec = [
    #     ### 
    #     # "61-64", # LA to SF, 348 miles, b=50, t=200 seems ok.
    #     # "61-41", # LA to Phoenix 356 miles, b=100, t=200 seems ok.

    #     "121-122", # Jacksonville to Miami, 328 miles

    #     ## > 500 mile OD network
    #     "61-491",  # From LA to Salt Lake City ## b=100, t=200 failed, 
    #     "131-472", # From Atlanta to Nashville 579 miles
    #     "484-485", # From Dallas to El Paso, 554 miles, b=100, t=200 seems ok.
    #     "261-452", # From Detroit to Greenville, 519 miles
    #     "486-485", # From Houston to El paso, 659 miles


    #     "61-531", # From LA to Seattle, 959 miles

    #     
    # ]
    return od_region_vec
end


function get_alg_label(alg_name)
    alg_label_dict = Dict(zip(g_alg_names, g_alg_xlabels))
    return alg_label_dict[alg_name]
end

function set_result_continent!(sym)
    @assert sym in [:us, :eu]
    global g_result_continent  = sym
    if g_result_continent == :us
        global g_cfo_result_path = joinpath(ep.k_data_path, "results", "cfo-results.oldmap.2025-05-15.jld2") 
    elseif  g_result_continent == :eu
        global g_cfo_result_path = joinpath(ep.k_data_path, "results", "cfo-results.eu.2025-05-15.jld2") 
    end
end
    


