"""
The datasets for EIA. This is the initial dataset that we chose.
"""

const g_eia_carbon_region_dict = Dict(
	"CAL"=>"California",
	"CAR"=>"Carolinas",
	"CENT"=>"Central",
	"FLA"=>"Florida",
	"MIDA"=>"Mid-Atlantic",
	"MIDW"=>"Midwest",
	"NE"=>"New England",
	"NY"=>"New York",
	"NW"=>"Northwest",
	"SE"=>"Southeast",
	"SW"=> "Southwest",
	"TEN"=>"Tennessee",
	"TEX"=>"Texas",
)
const g_eia_carbon_regions = collect(keys(g_eia_carbon_region_dict))

const g_energy_type_dict = Dict(
	"WAT"=>"Hydro",
	"SUN"=>"Solar",
	"WND"=>"Wind",
	"COL"=>"Coal",
	"NG" =>"NatureGas",
	"NUC"=>"Nuclear",
	"OTH"=>"Other",
	"OIL"=>"Petroleum",
	"ALL"=>"Total",
)

const g_renewable_types = ["WAT", "SUN", "WND"]
const g_pollution_types = ["COL", "NG", "OIL"]
const g_all_energy_types = ["COL", "NG", "OIL", "WAT", "SUN", "WND", "NUC", "OTH"]

# Refer to https=>//www.eia.gov/tools/glossary/index.php 
# Refer to https=>//www.eia.gov/electricity/gridmonitor/dashboard/electric_overview/US48/US48
# Refer to https=>//www.eia.gov/state/
const g_eia_region2state_dict = Dict(
	"CAL"=>["California"],
	"FLA"=>["Florida"],
	"NY"=>["New York"],
	"TEX"=>["Texas"],
	"TEN"=>["Tennessee"],
	"SE"=>[ "Georgia", "Alabama"],
	"SW"=> ["Arizona", "New Mexico"],
	"CAR"=>["North Carolina", "South Carolina"],
	"NE"=>["Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont"],

	"NW"=>["Oregon", "Washington", "Idaho", "Montana", "Wyoming", "Colorado", "Nevada", "Utah"],

	"CENT"=>["Oklahoma", "Kansas","Nebraska","North Dakota", "South Dakota",],  #, ],

	"MIDA"=>["Delaware", "District of Columbia", "Maryland", "New Jersey", "Pennsylvania", "Virginia", "West Virginia", "Ohio", "Kentucky"],

	"MIDW"=>["Illinois", "Louisiana", "Arkansas","Indiana", "Iowa",  "Michigan", "Minnesota", "Missouri",   "Wisconsin", "Mississippi"],
)

global g_eia_state2region_dict::OrderedDict{String, String} = Dict{String, String}()

function init_cs_eia()
    for (reg, states) in g_eia_region2state_dict
        for state in states
            global g_eia_state2region_dict[state] = reg
        end
    end
end

get_carbon_dict_predict(dataset::EiaDataset) = get_carbon_dict(dataset)

function get_carbon_dict(dataset::EiaDataset)
    carbon_regions = g_eia_carbon_regions
    carbon_dict = OrderedDict(carbon_regions .=> eia_get_carbon_vec.(carbon_regions))
    carbon_dict1 = OrderedDict(
        "MidCase" => carbon_dict,
    )
    return carbon_dict1
end

function latlon2carbonregion(lat::Real, lon::Real, region::Symbol, carbon_dataset::EiaDataset, region_box_dict)
    @assert (region == :us)

    state = latlon2subregion(lat, lon, Symbol(region))

    if !(state in g_statesUS48)
        return ""
    end

    reg = g_eia_state2region_dict[state]

    return reg
    
end

#####################
# carbon functions for EIA dataset


function eia_get_carbon_df(region::String)
    datapath = joinpath(g_carbon_data_dir, "$region.csv")
    df = DataFrame(CSV.File(datapath; dateformat="yyyy-mm-dd HH:MM:SS"))
    filter!(row -> row.timestamp >= DateTime(2018, 7, 1), df)
    df[!,:ALL] = convert.(Float64, df[!,:ALL])
    # filter!(row -> row.timestamp >= DateTime(2020, 1, 1), df)
    return df
end


"""
Get the price of carbon for each region.

renewable_mul: the multiplier for renewable energy, we may change it because we may vary the renewable percentage in the grid.
"""
function eia_get_carbon_vec(region::String, renewable_mul::Float64 = 1.0)
    df = eia_get_carbon_df(region)
    n_data = length(df.ALL)
    price_vec = zeros(Float64, n_data)
    for type in g_renewable_types
        data_vec = getproperty(df, Symbol(type))
        df[!, Symbol(type)] = data_vec * renewable_mul
    end

    for row in eachrow(df)
        all_val = 0.0
        for type in g_all_energy_types
            val = row[Symbol(type)]
            if !ismissing(val)
                all_val += val
            end
        end
        row.ALL = all_val
    end

    for type in g_pollution_types
        price = g_unit_co2_dict[type]
        data_vec = getproperty(df, Symbol(type))
        data_vec = data_vec * price
        price_vec += [ismissing(x) ? 0.0 : x for x in data_vec]
    end
    # format = DateFormat("yyyy-mm-dd HH:MM:SS")
    # timestamp0 = DateTime.(df.timestamp, (format,)) 
    timestamp0 = df.timestamp
    # convert UTC to UTC-8 (Pacific time zone)
    timestamp = timestamp0 .- Hour(8)
    price_vec = price_vec ./ df.ALL
    # data = (datetime=timestamp, price=price_vec)
    # ta0 = TS.TimeArray(data; timestamp = :datetime)

    ## Note that the carbon data may miss some data, so we need to fill in those data.
    full_ts_vec = collect(first(timestamp):Hour(1):last(timestamp))
    full_price_vec = fill(NaN, length(full_ts_vec))
    for ts in timestamp
        full_idx = searchsortedfirst(full_ts_vec, ts)
        idx = searchsortedfirst(timestamp, ts)
        full_price_vec[full_idx] = price_vec[idx]
    end
    full_price_vec = collect(full_price_vec)

    # ta = TS.TimeArray((datetime=full_ts_vec, price=full_price_vec); timestamp = :datetime)
    ta = TS.TimeArray(full_ts_vec, full_price_vec)

    return ta
end
# 

function eia_get_renewable_penetration()
    renewable_gen = 0.0
    total_gen = 0.0
    for region in g_eia_carbon_regions
        df = eia_get_carbon_df(region)
        for type in g_renewable_types
            for x in df[!, Symbol(type)]
                if ismissing(x)
                    continue
                end
                renewable_gen += x
            end
        end
        total_gen += sum(df.ALL)
    end
    return renewable_gen / total_gen
end
 
function eia_renewable_multiplier_to_ratio(multiplier::Real)
    current_ratio = eia_get_renewable_penetration()
    ratio = (multiplier * current_ratio) / (multiplier * current_ratio + 1 - current_ratio)
    return ratio
end
 
"""
Compute the renewable multiplier given the target percentage of renewable energy.
"""
function eia_renewable_ratio_to_multiplier(ratio::Real)
    
    current_ratio = eia_get_renewable_penetration()
    target_ratio = ratio
    # @show renewable_gen total_gen 
    tc = current_ratio * target_ratio
    multiplier = (target_ratio - tc) / (current_ratio - tc)
 
    return multiplier
end