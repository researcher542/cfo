mutable struct SpeedPlanPoint <: AbstractNode
	lat::Float64	
	lon::Float64
	ele::Float64
	speed::Float64
	distance::Float64
end

SpeedPlanPoint(node::AbstractNode, speed::Float64, distance::Float64) = SpeedPlanPoint(node.lat, node.lon, node.ele, speed, distance)

function Base.getproperty(pt::SpeedPlanPoint, name::Symbol)
    if name == :rpm
        return speed2rpm(pt.speed)
    elseif name == :duration
        return pt.speed != 0.0 ? pt.distance / pt.speed : 0.0
    else
        return getfield(pt, name)
    end
end

mutable struct Step 
	summary::Summary	
	speedplan::Vector{SpeedPlanPoint}
	polyline::String
end

function Base.:(==)(s1::Step, s2::Step)
    return (s1.polyline == s2.polyline) && (s1.summary == s2.summary)
end

Base.show(io::IO, s::Step) = show(io, s.summary)

Step() = Step(Summary(), [], "")

"""
The total duration of a plan
"""
planduration(plan::SpeedPlanPoint)::Float64 = plan.speed != 0.0 ? plan.distance / plan.speed : 0.0
planduration(plans::Vector{SpeedPlanPoint}) = sum(planduration.(plans))
planduration(step::Step) = planduration(step.speedplan)

"""
check if a step is highway, or should it be optimized
"""
function ishighway(step::Step)
	dis = step.summary.distance
	time = step.summary.duration
	spd = dis / time
	flag1 = dis > g_min_highway_distance
	flag2 = spd > g_min_highway_speed
	return flag1 && flag2
end