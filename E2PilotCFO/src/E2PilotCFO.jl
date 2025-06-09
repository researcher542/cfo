module E2PilotCFO


using Parameters, Dates, Graphs, Infiltrator, SparseArrays, TimerOutputs, Printf, JuMP, DataStructures, StaticArrays, Roots, CSV, DataFrames, Statistics, Geodesy, EzXML, JLD2, Distances, ProgressMeter, FastPriorityQueues, Interpolations, AngleBetweenVectors 

import OffsetArrays as OA
import TimeSeries as TS
import NearestNeighbors as NN
import LinearAlgebra as LA
import BlackBoxOptim as BBO

## Import the core functions
import E2PilotCore as ep
import E2PilotCore: k_data_path, g_diesel_j_per_liter, g_diesel_j_per_liter, g_liter_per_gallon, g_j_per_kwh, g_net_data_path
import E2PilotCore: AbstractObjective, ObjEnergy, ObjCarbon, ObjTime
import E2PilotCore: Node, Way, Network, AbsNet
import E2PilotCore: AbstractVehicle, PWL
import E2PilotCore: outneighbors, inneighbors, all_neighbors, get_minmax_speed, getway, getnode, output_power, g_to, energy_cost_on_road, closest_node, latlon2subregion, is_endpoint, distance3d, j2kwh, j2liter, set_task_tid, max_lambda, pos, e2map, ArrayDict
import E2PilotCore: Step, Astar, shortest_path, multi_paso, paso, fastest_path, minmax_t, get_minmax_t, augmentedtime
import E2PilotCore: AbstractPredictionMode, PredictNoise, PredictML, PredictPerfect
import E2PilotCore: FAF, ETIS, AbstractODSet


## Define some types
abstract type AbstractCarbonDataset end
struct CambiumDataset <: AbstractCarbonDataset end
struct CarbonCastDataset <: AbstractCarbonDataset end
struct EiaDataset <: AbstractCarbonDataset end
struct ElecMapDataset <: AbstractCarbonDataset end ## data set from ElecMapDataset


include("config.jl")
include("ev.jl")
include("prob.jl")
include("charge/charge.jl")
include("cs_net/cs_net.jl")
include("algs/algs.jl")
include("utils/utils.jl")
# include("approx/approx.jl")
include("results/results.jl")
include("testdata/testdata.jl")

global g_project_name::Symbol = :CFO

function init()
    @info "Initilizing E2PilotCFO" 
    ep.init()
    init_cs_eia()
end

end # module E2PilotCFO
