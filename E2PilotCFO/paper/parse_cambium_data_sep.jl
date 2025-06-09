"""
Parse the cambium data from the separate csv file and save it to a JLD2 file.
"""

import E2PilotCFO as cfo
import XLSX
import TimeSeries as TS
using Dates
using Infiltrator
using DataStructures
using JLD2
using CSV, DataFrames

sf = safehouse


datapath = joinpath(cfo.g_carbon_data_dir, "cambium", "hourly_state")

function get_file_path(state, scenario, year)
    filename = "Cambium22_$(scenario)_hourly_$(state)_$(year).csv"
    return joinpath(datapath, filename)
end

function get_ts_val_vec(state, scenario, year)
    file_path = get_file_path(state, scenario, year)
    df = DataFrame(CSV.File(file_path; header=6))
    ts_vec0 = df[!, "timestamp"]
    fmt = DateFormat("yyyy-mm-dd HH:MM:SS")
    ts_vec = [DateTime(x, fmt) for x in ts_vec0]
    value_vec = df[!, value_name]
    ts_val_vec = [(ts, val) for (ts, val) in zip(ts_vec, value_vec)]
    return ts_val_vec
end



value_name = ["aer_load_co2e", "aer_gen_co2e", "srmer_co2e", "lrmer_co2e"][4]
scenario_vec = ["HighRECost", "Electrification", "MidCase", "MidCase100by2035", "MidCase95by2050"]
state_vec = cfo.g_us_abbr48_vec

function parse_cambium_data()
    data_dict = OrderedDict()
    for scenario in scenario_vec
        data_dict[scenario] = OrderedDict()
        for state in state_vec
            @show scenario, state
            ts_val_vec = []
            for year in [2024, 2026, 2028, 2030, 2035, 2040, 2045, 2050]
                ts_val_vec1 = get_ts_val_vec(state, scenario, year)
                ts_val_vec = vcat(ts_val_vec, ts_val_vec1)
            end
            sort!(ts_val_vec; by=x->x[1])
            ts_vec = [x[1] for x in ts_val_vec]
            value_vec = [ (x[2] + 1e-6)  for x in ts_val_vec]
            value_vec = value_vec ./ 1000.0
            data_dict[scenario][state] = TS.TimeArray(ts_vec, value_vec)
        end
    end

    data_dict2 = OrderedDict(
        scenario => OrderedDict(
               state => data_dict[scenario][state]
                   for state in keys(data_dict[scenario])
           )
           for scenario in keys(data_dict)
    )

    return data_dict2
end

data_dict2 = parse_cambium_data()
carbon_dict_path = cfo.get_carbon_dict_path(value_name)
@time JLD2.save_object(carbon_dict_path, data_dict2)

# ta = TS.TimeArray(ts_vec, value_vec)

#function parse_combium_data(scenario, )
# end

# @time data_dict = parse_combium_data(xf)
# @time JLD2.save_object(cfo.g_carbon_dict_path, data_dict)

