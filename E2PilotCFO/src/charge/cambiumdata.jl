

get_region_box_dict(::AbstractCarbonDataset) = OrderedDict()

function get_carbon_dict_predict(dataset::CambiumDataset, noise_ratio::Float64 = 0.1)
    carbon_dict = get_carbon_dict(dataset) 
    for scenario in keys(carbon_dict)
        for reg in keys(carbon_dict[scenario])
            ta = carbon_dict[scenario][reg]
            n_data = length(TS.timestamp(ta))
            noise = randn(n_data) * noise_ratio .* TS.values(ta)
            new_ta = TS.TimeArray(TS.timestamp(ta), TS.values(ta) .+ noise)
            carbon_dict[scenario][reg] = new_ta
        end
    end 
    return carbon_dict
end

global g_carbon_dict_suffix_vec = ["aer_load_co2e", "aer_gen_co2e", "srmer_co2e"]
global g_default_carbon_dict_suffix = g_carbon_dict_suffix_vec[1]
function get_carbon_dict(dataset::CambiumDataset, suffix::String = g_default_carbon_dict_suffix)
    carbon_dict_path = get_carbon_dict_path(suffix)
    if !isfile(carbon_dict_path)
        @warn "Carbon data not found, please run the script 'parse_combium_data.jl' to parse the data from Combium."
        return missing
    end
    carbon_dict = JLD2.load_object(carbon_dict_path)
    # carbon_dict = OrderedDict(g_eia_carbon_regions .=> get_carbon_vec.(g_eia_carbon_regions, renewable_mul))
    return carbon_dict
end

function latlon2carbonregion(lat::Real, lon::Real, region::Symbol, carbon_dataset::CambiumDataset, region_box_dict)
    state = latlon2subregion(lat, lon, Symbol(region))
    if !(state in g_statesUS48)
        return ""
    end
    reg = g_us_state_abbr_dict[state]
    return reg
end