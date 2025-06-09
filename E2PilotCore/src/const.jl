"""
define some global constants
"""


const g_to = TimerOutput()

const g_theta_range = (-0.1:5e-4:0.1) # unitless
const g_spd_range = (1.0:0.1:35.0) # in m/s
const g_min_theta = -0.1
const g_max_theta = 0.1
const g_min_spd = 1.0
const g_max_spd = 35.0

# The minimium require to optimize, we will skip the shorter roads to make our assumption valid.
const g_min_highway_distance = 5_000.0 
const g_min_highway_speed =  (48/3.6) # 48 kmh
const g_max_highway_speed = (120/3.6)

# The minimium distance for two nodes. Smaller distance will cause large grade due to the resolution of the elevation data, while larger distance will skip the up-and-down roads and make our algorithm inaccurate.
const k_min_nodes_distance = 2_000.0 

const k_root_dir = normpath(joinpath(@__DIR__, "../.."))
const k_data_path = normpath(joinpath(k_root_dir), "data")
const k_fig_dir = normpath(joinpath(k_root_dir), "figs")

