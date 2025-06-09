

get_carbon_dict_predict(dataset::ElecMapDataset) = get_carbon_dict(dataset)

function get_region_vec(dataset::ElecMapDataset)
    reg_vec = ["AL", "AT", "BA", "BE", "BG", "BY", "CH", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GB", "GR", "HR", "HU", "IE", "IS", "IT", "LT", "LU", "LV", "MD", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK"]
    return reg_vec
end

function get_fullname2region_dict(dataset::ElecMapDataset)
    region_dict = get_region2fullname_dict(dataset)
    abbr_dict = OrderedDict(
        val => key for (key, val) in region_dict
    )
    abbr_dict["Czech Republic"] = "CZ"
    abbr_dict["The former Yugoslav Republic of Macedonia"] = "MK"
    return abbr_dict 
end

"""
Map the abbr to region name
"""
function get_region2fullname_dict(dataset::ElecMapDataset)
    region_dict = Dict(
        # Western Europe
        "AT" => "Austria",
        "BE" => "Belgium",
        "FR" => "France",
        "DE" => "Germany",
        "LI" => "Liechtenstein",
        "LU" => "Luxembourg",
        "MC" => "Monaco",
        "NL" => "Netherlands",
        "CH" => "Switzerland",
    
        # Northern Europe
        "DK" => "Denmark",
        "EE" => "Estonia",
        "FI" => "Finland",
        "IS" => "Iceland",
        "IE" => "Ireland",
        "LV" => "Latvia",
        "LT" => "Lithuania",
        "NO" => "Norway",
        "SE" => "Sweden",
        # "GB" => "United Kingdom",
        "GB" => "Great Britain",
    
        # Southern Europe
        "AL" => "Albania",
        "AD" => "Andorra",
        "BA" => "Bosnia and Herzegovina",
        "HR" => "Croatia",
        "CY" => "Cyprus",
        "GR" => "Greece",
        "IT" => "Italy",
        "MT" => "Malta",
        "ME" => "Montenegro",
        "MK" => "North Macedonia",
        "PT" => "Portugal",
        "SM" => "San Marino",
        "RS" => "Serbia",
        "SI" => "Slovenia",
        "ES" => "Spain",
        "VA" => "Vatican City",
    
        # Eastern Europe
        "BY" => "Belarus",
        "BG" => "Bulgaria",
        # "CZ" => "Czech Republic",
        "CZ" => "Czechia",
        "HU" => "Hungary",
        "MD" => "Moldova",
        "PL" => "Poland",
        "RO" => "Romania",
        "RU" => "Russia",
        "SK" => "Slovakia",
        "UA" => "Ukraine",
    
        # Transcontinental (Europe/Asia)
        "AM" => "Armenia",
        "AZ" => "Azerbaijan",
        "GE" => "Georgia",
        "KZ" => "Kazakhstan",
        "TR" => "Turkey"
    ) 
    return region_dict
end

# function get_region_box_dict(dataset::ElecMapDataset)
#    reg_dict = OrderedDict(
#     )
#     
# end

function get_carbon_data_path(dataset::ElecMapDataset, region::String)
    # datapath = joinpath(g_carbon_data_dir, "elecmap", "$region.csv")
    folder = joinpath(g_carbon_data_dir, "elecmap")
    path = joinpath(folder, "$(region)_2024_hourly.csv")
    if !isfile(path)
        @warn("The file $path does not exist.")
        return missing
    end
    return path
end

"""
"""
function get_carbon_data(dataset::ElecMapDataset, region::String)
    path = get_carbon_data_path(dataset, region)
    df = DataFrame(CSV.File(path))
    zone_id = df[1, "Zone id"]
    @assert zone_id == region

    check_flag = true
    if check_flag
        country = df[1, "Country"] 
        region_dict = get_region2fullname_dict(dataset)
        if country != region_dict[region]
            @warn("The country $country does not match the region $region. In dict: $(region_dict[region])")
        end
    end
    # df1 = unique(sort(df, :datetime), :datetime)
    # key = "Carbon intensity gCO₂eq/kWh (direct)"
    key = "Carbon intensity gCO₂eq/kWh (Life cycle)"
    ci_vec = df[!, key] ./ 1000.0
    ts_vec = [DateTime(ts, "yyyy-mm-dd HH:MM:SS") for ts in df[!, "Datetime (UTC)"]]
    ta = TS.TimeArray(ts_vec, ci_vec)
    return ta
end


function get_carbon_dict(dataset::ElecMapDataset)
    region_vec = get_region_vec(dataset)
    d = OrderedDict(
        region_vec .=> [get_carbon_data(dataset, region) for region in region_vec]
    )

    return OrderedDict(
        "MidCase" => d
    )

end

function latlon2carbonregion(lat::Real, lon::Real, region::Symbol, carbon_dataset::ElecMapDataset, region_box_dict)
    # @show region
    @assert (region == :eu)
    fullname = latlon2subregion(lat, lon, region, 0)
    if isnothing(fullname)
        return ""
    end
    dict = get_fullname2region_dict(carbon_dataset)
    reg = dict[fullname]
    if !(reg in get_region_vec(carbon_dataset))
        # @warn("The region $reg is not in the region vector.")
        return ""
    end
    return reg
end