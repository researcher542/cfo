"""
"""


# To resolve the recursive definition 
abstract type AbstractVehicle end

abstract type AbstractObjective end
struct ObjEnergy <: AbstractObjective end
struct ObjTime <: AbstractObjective end
struct ObjEmission <: AbstractObjective end
struct ObjCarbon <: AbstractObjective end

## different prediction modes. Currently used 
abstract type AbstractPredictionMode end
struct PredictNoise <: AbstractPredictionMode end #add noise to the prediction
struct PredictML <: AbstractPredictionMode end #use machine learning to predict
struct PredictPerfect <: AbstractPredictionMode end #perfect prediction, directly use the true value


include("solgrid.jl")
include("vehicle/vehicle.jl")
include("speedplan.jl")
include("route/routeplan.jl")
include("cost_sum.jl")
include("paso.jl")
include("multi_paso.jl")
include("augmenttime_raw.jl")
include("optimal_speed.jl")

#module E2Plan
#end