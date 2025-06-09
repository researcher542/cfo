"""
Emission related functions. 
Used for journal extension of the BuildSys'21
"""


"""
The parameters for emission functions (k1, k2, offset) tuple
"""
const k_emission_param_vec = [
    #(6128.0, 19500.0, 5.896e-16),
    #(2807.0, 16821.0, 5.896e-16),
    #(10710.0, 20690.0, 7.389e-16),
    (2703.0, 16550.0, 5.896e-16),
    (12096.0, 21266.0, 2.0), # const offset to make sure multi-injection is always better.
]
g_switching_speed_vec = [3719.0,]

const k_min_rpm = 500.0
const k_max_rpm = 6500.0


"""
paramter: 
    engine_speed: in RPM

return:
    emission rate: in ppm/s [??Not sure what is the unit]
"""
function emission_function(engine_speed::Real, grade::Real, k1::Float64, k2::Float64, offset::Float64)
    coeff = grade2emission_coeef(grade)
    val = k1*exp(-k2/engine_speed) + offset
    return coeff*val
end

"""
Find the speed range [s_{i-1}, s_i] and corresponding parameters and the then compute the emission rate
"""
function emission_function(engine_speed::Real, grade::Real, param_vec::Vector, swithing_speed_vec::Vector)
    idx = searchsortedfirst(swithing_speed_vec, engine_speed)
    param = param_vec[idx]
    return emission_function(engine_speed, grade, param...)
end


"""
The derivative of the emission function
"""
function emission_function_d(engine_speed::Float64, grade::Float64, k1::Float64, k2::Float64)
    coeff = grade2emission_coeef(grade)
    val = k1*k2*exp(-k2/engine_speed)/(engine_speed^2)
    return coeff*val
end

"""
Emission function also depends on grade. 
So we need to adjust the emission function by a coefficient.
But I don't know how it works, just copied here.
"""
function grade2emission_coeef(grade::Float64)::Float64
    #return 1.0
    if grade >= 0 
        return 1+70*grade
    elseif grade >= -0.024
        return 1+35*grade
    else
        return 0.16
    end
end

"""
Find the get generailized right point p of a function:
g(x) = k1 exp^{-k2/x} + offset. given (x1, y1) below g(x)
the point is given by the equation: g'(p) = (g(p) - y1) / (p-x1)
---
parameters:

    grade, k1, k2, offset: used for emission function
    x1,y1: the point below g(x)
"""
function get_generailized_right_point(grade::Float64, k1::Float64, k2::Float64, offset::Float64, x1::Float64, y1::Float64)
    # The function to find zero
    function func(p)
        g_d = emission_function_d(p, grade, k1, k2)
        g_p = emission_function(p, grade, k1, k2, offset)
        rhs = (g_p - y1) / (p - x1)
        return g_d - rhs
    end 

    p0 = x1+1e0
    p = find_zero(func, p0; verbose=false)
    return p
end

"""
Get the objective and solution of the convex optimization problem in the paper Theorem~6 equation (12)

engine_speed: average speed

return:
---
    objective: the optimal objective 
    β: vector: β_i is the fraction of time the speed_i is used.
    spd_vec
"""
function emission_objective(grade::Float64, engine_speed::Float64, param_vec::Vector, switching_speed_vec::Vector)
    n_piece = length(param_vec)
    beta_vec = zeros(Float64, n_piece)
    speed_vec = zeros(Float64, n_piece)
    if isempty(switching_speed_vec) || (engine_speed <= switching_speed_vec[1])
       beta_vec[1] = 1 
       speed_vec[1] = engine_speed
    else
        # idx must be >= 1 
        idx = searchsortedfirst(switching_speed_vec, engine_speed) - 1
        s_im1 = switching_speed_vec[idx] 
        y1 = emission_function(s_im1, grade, param_vec[idx]...)
        p = get_generailized_right_point(grade, param_vec[idx+1]..., s_im1, y1)
        #@debug "p=$p, ω=$engine_speed"
        if p > engine_speed
            si = idx == n_piece-1 ? k_max_rpm : switching_speed_vec[idx+1]
            ri = min(p, si)
            # @debug "ri=$ri"
            beta1 = (ri-engine_speed)/(ri-s_im1)
            beta2 = 1 - beta1
            beta_vec[idx] = beta1
            beta_vec[idx+1] = beta2
            speed_vec[idx] = s_im1
            speed_vec[idx+1] = ri
            # @debug "two points in emission"
        else
            speed_vec[idx+1] = engine_speed
            beta_vec[idx+1] = 1
        end
    end
    # 
    obj = 0.0
    for i in 1:n_piece
        if beta_vec[i] == 0 continue end
        betai = beta_vec[i]
        spd = speed_vec[i]
        param = param_vec[i]
        val = emission_function(spd, grade, param...)
        obj += betai * val
    end
    #@debug "spd:$engine_speed, beta: $beta_vec, obj: $obj"

    val1 = emission_function(engine_speed, grade, param_vec, switching_speed_vec)
    diff = val1 - obj
    avg_spd = sum(beta_vec .* speed_vec)
    diff2 = (avg_spd - engine_speed)
    if !(abs(diff2) <= 1e-10)
        @warn "average speed not accurate with diff: $diff2"
    end
    if !(diff ≈ 0)
        #@debug "For emission objective, rpm:$engine_speed diff:$diff < 0, calculated obj $obj, val1:$val1"
    end

    return obj, beta_vec, speed_vec
end

const k_rpm_per_mph = 80.0
const k_transmission_ratio = 4.0
const k_wheel_diam = 0.5
#k_rpm_per_ms = 1/(k_wheel_diam*π/60/k_transmission_ratio)
const k_rpm_per_ms = k_rpm_per_mph * g_mph_per_ms
speed2rpm(speed::Real) = (speed*k_rpm_per_ms)
ms2rpm = speed2rpm
rpm2speed(rpm::Real) = (rpm/k_rpm_per_ms)
rpm2ms = rpm2speed

"""
Given theta and speed, compute the emission objective 
"""
function objective(veh::AbstractVehicle, theta::Real, speed::Real, objtype::ObjEmission; param_vec=k_emission_param_vec, swi_spd_vec=g_switching_speed_vec)
    engine_speed = speed2rpm(speed)
    val = emission_objective(theta, engine_speed, param_vec, swi_spd_vec)[1]
    #@debug "objective called with spd=$swi_spd_vec"
    return val
end


emission_cost_v(veh::Vehicle, speedplan::Vector{SpeedPlanPoint}, param_vec::Vector, swiching_speed_vec::Vector) = emission_cost_v(veh, speedplan, [pt.speed for pt in speedplan[1:end-1]] , param_vec, swiching_speed_vec)
emission_cost_v(veh::Vehicle, step::Step, param_vec::Vector, swiching_speed_vec::Vector) = emission_cost_v(veh, step.speedplan , param_vec, swiching_speed_vec)

function emission_cost_v(veh::Vehicle, nodes::Vector{S}, spd_v::Vector, param_vec::Vector, swiching_speed_vec::Vector) where S <: AbstractNode
    dis_v = distance3d.(nodes[1:end-1], nodes[2:end])
	grade_v = grades(nodes)
	div(a,b) = (b == 0) ? 0.0 : a/b
	time_v = div.(dis_v, spd_v)
    rpm_v = speed2rpm.(spd_v)
    # should use emission_function instead of emission objective:
    # In emission objective, we will seprate the speed. If we call it twice, then we will get sub-optimal solution.
    emisstion_rate_v = emission_function.(rpm_v, grade_v, (param_vec,), (swiching_speed_vec,))
    #emisstion_rate_v = emission_objective.(grade_v, rpm_v, (param_vec,), (swiching_speed_vec,))
    cost_v = [time_v[i] * emisstion_rate_v[i] for i in 1:length(time_v) ]
    return cost_v
end

total_emission_cost(veh::Vehicle, speedplan::Vector{SpeedPlanPoint}, param_vec::Vector, swiching_speed_vec::Vector) = total_emission_cost(veh, speedplan, [pt.speed for pt in speedplan[1:end-1]] , param_vec, swiching_speed_vec)
total_emission_cost(veh::Vehicle, step::Step, param_vec::Vector, swiching_speed_vec::Vector) = total_emission_cost(veh, step.speedplan , param_vec, swiching_speed_vec)

function total_emission_cost(veh::Vehicle, nodes::Vector{S}, spd_v::Vector, param_vec::Vector, swiching_speed_vec::Vector) where S <: AbstractNode
    cost_v = emission_cost_v(veh, nodes, spd_v, param_vec, swiching_speed_vec) 
	cost = sum(cost_v)
    return cost
end

function step2driving_cycle(step::Step)
    speed_v = Vector{Vector{Float64}}()
    for pt in step.speedplan
        # duration = pt.speed != 0 ? pt.distance / pt.speed : 0
        n_seconds = Int64(round(pt.duration))
        rpm = pt.rpm
        push!(speed_v, rpm*ones(n_seconds))
    end
    return vcat(speed_v...)
end

function step2grade_cycle(step::Step)
    grade_v_v = Vector{Vector{Float64}}()
    plan = step.speedplan
    grade_v = grades(plan)
    for i in 1:length(plan)-1
        pt = plan[i]
        n_seconds = Int64(round(pt.duration))
        θ = grade_v[i]
        push!(grade_v_v, θ*ones(n_seconds))
    end
    pt = plan[end]
    n_seconds = Int64(round(pt.duration))
    push!(grade_v_v, grade_v[end]*ones(n_seconds))

    return vcat(grade_v_v...)
end

"""get the emission cost from driving cycle, one second per line."""
function emisstion_cost_davr(step::Step, index_vec::Vector; param_vec::Vector, kwargs...)
    params = [param_vec[idx+1] for idx in index_vec]
    dc_rpm_vec = step2driving_cycle(step)
    dc_grade_vec = step2grade_cycle(step)
    
    emission_cost_v = [
        emission_function(dc_rpm_vec[i], dc_grade_vec[i], params[i]...) for i in 1:length(dc_rpm_vec)
    ] 
    return sum(emission_cost_v)
end