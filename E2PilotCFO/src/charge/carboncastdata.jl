"""
Use the data from CarbonCast: https://github.com/carbonfirst/CarbonCast with branch `v3.0`     
"""



function get_region_vec(dataset::CarbonCastDataset)
    return ["AECI", "AZPS", "BPAT", "CISO", "DUK", "EPE", "ERCO", "FPL", "ISNE", "LDWP", "MISO", "NEVP",
        "NWMT", "NYIS", "PACE", "PJM", "SC", "SCEG", "SOCO", "SRP",
        "TIDC", "TVA", "WACM"]
end

"""
"""
function get_region_box_dict(dataset::CarbonCastDataset)
    ## copied from CarbonCast/src/weather/ds084.1_control.ctl
    bounding_box_str = """
    CISO: nlat=42 slat=32 wlon=-124.75 elon=-113.5
    PJM: nlat=43 slat=34.25 wlon=-91 elon=-73.5
    ERCOT: nlat=36.5 slat=25.25 wlon=-104.5 elon=-93.25
    ISNE: nlat=48 slat=40 wlon=-74.25 elon=-66.5
    MISO: nlat=50.00 slat=28.50 wlon=-107.75 elon=-81.75
    BPAT: nlat=49.50 slat=39.50 wlon=-125.25 elon=-105.50
    SWPP: nlat=49.50 slat=30.25 wlon=-107.75 elon=-89.50
    SOCO: nlat=35.50 slat=29.25 wlon=-90.50 elon=-80.25 
    FPL: nlat=31.25 slat=24.00 wlon=-83.50 elon=-79.50 
    NYISO: nlat=45.50 slat=40.00 wlon=-80.25 elon=-71.25
    BANC: nlat=41.75 slat=37.00 wlon=-124.00 elon=-120.00
    LDWP: nlat=38.00 slat=33.25 wlon=-119.00 elon=-117.00
    TIDC: nlat=38.25 slat=36.75 wlon=-121.75 elon=-119.75
    DUK: nlat=37.00 slat=33.00 wlon=-84.75 elon=-77.75 
    SC: nlat=35.25 slat=31.50 wlon=-82.75 elon=-78.00 
    SCEG: nlat=35.25 slat=31.50 wlon=-83.00 elon=-78.75 
    SPA: nlat=40.75 slat=34.25 wlon=-98.00 elon=-89.00 
    FMPP: nlat=30.75 slat=24.00 wlon=-83.00 elon=-79.50 
    FPC: nlat=31.25 slat=25.75 wlon=-86.50 elon=-80.00 
    TAL: nlat=31.25 slat=29.75 wlon=-84.75 elon=-83.50 
    TEC: nlat=29.00 slat=27.00 wlon=-83.25 elon=-81.25 
    AECI: nlat=41.75 slat=34.25 wlon=-98.50 elon=-88.50 
    LGEE: nlat=39.50 slat=36.00 wlon=-89.75 elon=-82.25 
    DOPD: nlat=49.50 slat=46.75 wlon=-120.75 elon=-118.25 
    GCPD: nlat=48.50 slat=46.25 wlon=-120.50 elon=-118.50 
    GRID: nlat=46.25 slat=44.75 wlon=-119.75 elon=-118.25 
    IPCO: nlat=47.25 slat=41.50 wlon=-120.50 elon=-111.00 
    NEVP: nlat=42.50 slat=34.50 wlon=-122.00 elon=-111.00 
    NWMT: nlat=49.50 slat=43.25 wlon=-116.50 elon=-103.50 
    PACE: nlat=45.50 slat=33.00 wlon=-115.75 elon=-104.25 
    PACW: nlat=47.50 slat=38.75 wlon=-124.75 elon=-115.75 
    PGE: nlat=46.50 slat=44.25 wlon=-124.25 elon=-121.25 
    PSCO: nlat=41.75 slat=35.75 wlon=-109.50 elon=-102.00 
    PSEI: nlat=49.50 slat=45.75 wlon=-123.75 elon=-119.75 
    SCL: nlat=48.25 slat=47.00 wlon=-123.00 elon=-121.75 
    TPWR: nlat=48.25 slat=45.75 wlon=-124.00 elon=-120.50 
    WACM: nlat=48.00  slat=35.50 wlon=-114.50 elon=-95.75 
    SOCO: nlat=35.50 slat=29.50 wlon=-90.50 elon=-80.25 
    AZPS: nlat=36.75 slat=30.75 wlon=-115.25 elon=-108.75 
    EPE: nlat=34.00  slat=26.75 wlon=-108.75 elon=-98.25 
    PNM: nlat=44.50 slat=30.75 wlon=-123.50 elon=-101.50 
    SRP: nlat=34.50 slat=32.00 wlon=-113.75 elon=-110.50 
    TEPC: nlat=36.75 slat=31.25 wlon=-115.25 elon=-110.00 
    WALC: nlat=44.00 slat=30.75 wlon=-124.25 elon=-105.00 
    TVA: nlat=38.00 slat=31.75 wlon=-90.75 elon=-81.25
    """

    reg_dict = OrderedDict()
    
    for line in split(bounding_box_str, "\n")
        if length(line) == 0
            continue
        end
        reg, latlon_str = split(line, ": ")

        latlon_vec = split(latlon_str, " ")
        filter!(x -> length(x) > 0, latlon_vec)
        # @show latlon_str latlon_vec 
        val_vec = (
            parse(Float64, split(latlon_vec[1], "=")[2]),
            parse(Float64, split(latlon_vec[2], "=")[2]),
            parse(Float64, split(latlon_vec[3], "=")[2]),
            parse(Float64, split(latlon_vec[4], "=")[2])
        )
        latlon_tuple = NamedTuple{(:nlat, :slat, :wlon, :elon)}(val_vec)
        reg_dict[reg] = latlon_tuple
    end
    reg_dict = OrderedDict(reg_dict...)
    region_vec = get_region_vec(dataset)
    for (key, val) in reg_dict
        if !(key in region_vec)
            delete!(reg_dict, key)
        end
    end
    return reg_dict
end


function get_carbon_data_path(dataset::CarbonCastDataset, region::String)
    folder = joinpath(g_carbon_data_dir, "carboncast", region)
    path = joinpath(folder, "$(region)_direct_96hr_CI_forecasts_0.csv")
    if !isfile(path)
        @warn("The file $path does not exist.")
        return missing
    end
    return path
end

"""
"""
function get_carbon_data(dataset::CarbonCastDataset, region::String, predict_flag::Bool)
    path = get_carbon_data_path(dataset, region)
    df = DataFrame(CSV.File(path))
    df1 = unique(sort(df, :datetime), :datetime)
    key = predict_flag ? "avg_carbon_intensity_forecast" : "carbon_intensity_actual"
    ci_vec = df1[!, key] ./ 1000.0
    ts_vec = df1[!, "datetime"]
    ta = TS.TimeArray(ts_vec, ci_vec)
    return ta
end

function _get_carbon_dict_carboncast(predict_flag::Bool)
    region_vec = get_region_vec(CarbonCastDataset())
    d = OrderedDict(
        region_vec .=> [get_carbon_data(CarbonCastDataset(), region, predict_flag) for region in region_vec]
    )

    return OrderedDict(
        "MidCase" => d
    )
end

function get_carbon_dict_predict(dataset::CarbonCastDataset)
    return _get_carbon_dict_carboncast(true)
end

function get_carbon_dict(dataset::CarbonCastDataset)
    return _get_carbon_dict_carboncast(false)
end

function latlon2carbonregion(lat::Real, lon::Real, region::Symbol, carbon_dataset::DT, region_box_dict) where DT <: Union{CarbonCastDataset,}
    @assert (region == :us || region == :eu)

    reg = ""
    for (key, val) in region_box_dict
        if val.nlat >= lat >= val.slat && val.wlon <= lon <= val.elon
            reg = key
            break
        end
    end

    function distance2center(val)
        lat_center = (val.nlat + val.slat) / 2
        lon_center = (val.wlon + val.elon) / 2
        latlon = LatLon(lat, lon)
        latlon_center = LatLon(lat_center, lon_center)
        return ep.distance2d(latlon, latlon_center)
    end

    ## if it is not in the dataset, then we might want to find the closest region
    if reg == ""
        value, reg = findmin(distance2center, region_box_dict)
        # @show reg value
    end

    return reg
end