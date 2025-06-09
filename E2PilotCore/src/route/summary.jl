import Base

mutable struct Summary
	cost::Float64
	duration::Float64
	distance::Float64
end

function Base.:(==)(s1::Summary, s2::Summary)
    return (s1.cost == s2.cost) && (s1.duration == s2.duration) && (s1.distance == s2.distance)
end


function Base.show(io::IO, s::Summary)
    cost = s.cost
    dis = s.distance/1000.0
    time = s.duration / 3600
    Base.print(io, "Summary: cost: $cost, distance: $dis km, time: $time hours")
end

function getspeed(s::Summary)::Float64
	return s.duration == 0 ? 0.0 : s.distance/s.duration
end

function Base.getproperty(s::Summary, name::Symbol)
	if name == :speed
		return getspeed(s)
	else
		return Base.getfield(s, name)
	end
end 

function Summary()
	return Summary(0, 0, 0)	
end

function Base.:+(sum1::Summary, sum2::Summary) 
	sum0 = Summary()
	for field in fieldnames(Summary)
		val1 = getfield(sum1, field)
		val2 = getfield(sum2, field)
		val =  val1 + val2
		setfield!(sum0, field, val)
	end
	return sum0
end
