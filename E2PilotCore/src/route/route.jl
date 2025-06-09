
include("summary.jl")
include("step.jl")

mutable struct Leg 
	steps::Vector{Step}
	summary::Summary	
end

mutable struct Route 
	legs::Vector{Leg}
	summary::Summary
	geometry::String  # The polyline of the total route
	name::String
	speedplan::Vector{SpeedPlanPoint}
	type::String
end

getsummary(road::Union{Route, Leg, Step}) = road.summary
getsteps(route::Route)::Vector{Step} = vcat(getsteps.(route.legs)...)
getsteps(leg::Leg)::Vector{Step} = leg.steps