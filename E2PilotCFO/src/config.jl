"""
"""

const g_max_wait_time = 3600.0 * 4.0  
const g_min_wait_time = 600.0
const g_min_charge_time = 0.0
const g_max_charge_time = 3600.0 * 2.0
g_cfo_proj_dir = joinpath(ep.k_root_dir, "E2PilotCFO")

const g_carbon_data_dir = joinpath(k_data_path, "carbon")
function get_carbon_dict_path(suffix::String)
    carbon_dict_path = joinpath(g_carbon_data_dir, "carbon_dict$suffix.jld2")
    return carbon_dict_path
end

const g_h2_kwh_per_kg = 33.3 # https://hypertextbook.com/facts/2005/MichelleFung.shtml
global g_electrolysis_eff = 0.80 
# the overall efficiency of electrolysis systems by converting water to H2 https://www.nrel.gov/docs/fy04osti/36705.pdf; Technology Brief: Analysis of Current-Day Commercial Electrolyzers 64%
# Hydrogen production by PEM water electrolysis – A review 80%
global g_fcev_eff = 0.6 # the efficiency of hydrogen fuel-cell EV A Detailed Vehicle Modeling & Simulation Study Quantifying Energy Consumption and Cost Reduction of Advanced Vehicle Technologies Through 2050; Fig.0-11

global g_hev_eff_improvement = 1.1 # the efficiency improvement of HEV compared to ICEV, https://vms.taps.anl.gov/analytics/md-hd-truck-future-technology-prediction/

global g_charge_eff = 0.9 # the efficiency of charging EVs.

# Refer to https://www.eia.gov/environment/emissions/co2_vol_mass.php
const g_kg_co2_per_gallon = 10.19
global g_kg_co2_per_j = g_kg_co2_per_gallon / (g_diesel_j_per_liter * g_liter_per_gallon)
global g_kg_co2_per_kwh = g_kg_co2_per_j * g_j_per_kwh

dieselenergy2carbon(diesel::Real) = (diesel*g_kg_co2_per_j) 

const g_default_scenario = "MidCase"
const g_default_scenario_key = "BAU"
const g_default_st_year = 2024

# CO2 emission per kwh, unit=> kg/kWh
# refer to https=>//www.eia.gov/tools/faqs/faq.php?id=74&t=11
const g_kg_per_pound = 0.453592
const g_unit_co2_dict = Dict(
	"WAT"=>0,
	"SUN"=>0,
	"WND"=>0,
	"COL"=>2.23*g_kg_per_pound,
	"NG" =>0.91*g_kg_per_pound,
	"NUC"=>0, # no emission during operation, but costly to mine.
	"OTH"=>0,
	"OIL"=> 2.13*g_kg_per_pound,
	"ALL"=> NaN,
)

g_global_warming_potential_dict = Dict(
    "CO2" => 1.0,
    "CH4" => 29.8,
    "N2O" => 273.0,
)

# global g_scenario_dict = Dict(
#     "BAU" => "Mid-case", 
#     "2050" => "95% Decarbonization by 2050", 
#     "2035" => "100% Decarbonization by 2035",
#     "elec" => "Electrification",
#     "highRenew" => "High Renewable Energy Costs",
# )

global g_scenario_dict = Dict(
    "BAU" => "MidCase", 
    "2050" => "MidCase95by2050", 
    "2035" => "MidCase100by2035",
    "elec" => "Electrification",
    "highRenew" => "HighRECost",
)

global g_selected_scenario_vec = ["highRenew", "BAU",  "2035"]

function get_grid_scenario_label(scenario)
    d = Dict(
        "BAU" => "Business as Usual", 
        "2050" => "95\\% Decarbonization by 2050", 
        "2035" => "100\\% Decarbonization by 2035"
    )
    return d[scenario]
end
    