



const g_default_rate = 1e6
# The piecewise function defined by a bunch of points

const g_perc_vec = SA[0.8, 0.85, 0.9, 0.95, 1]
const g_mul_vec = SA[1, 0.8, 0.6, 0.4, 0.2]

"""
β: the initial β

t: the charging time in seconds

B: the battery capacity in J

rate: the charging rate in W

return the final β after charging 
"""
function charge_function(β::Real, t::Real, B::Real, rate::Real = g_default_rate)
    n = length(g_perc_vec)
    if t <= 0.0
        return β
    end
    for i = 1:n
        perc = g_perc_vec[i]
        mul = g_mul_vec[i]
        required_time = (B*perc - β) / (mul*rate)
        if required_time < 0
            continue
        end
        if required_time < t
            t -= required_time
            β = perc*B
        else
            β += mul*rate * t
            break
        end
    end
    return β
end

"""
The inverse of the charge function, given final soc, return charging time
"""
function cf_inv(β0::Real, βf::Real, B::Real, rate::Real = g_default_rate)::Float64
    β = β0
    βf = clamp(βf, 0.0, B)
    n = length(g_perc_vec)
    t = 0.0
    for i = 1:n
        perc = g_perc_vec[i]
        mul = g_mul_vec[i]
        required_time = (B*perc - β) / (mul*rate)
        # 
        if required_time < 0
            continue
        end
        if B*perc < βf
            t += required_time
            β = perc*B
        else
            t += (βf - β) / (mul*rate)
            break
        end
    end
    return t
    # charged_β = βf - β
    # # TODO: need to be careful, not generic.
    # if βf < 0.8*B
    #     return charged_β/rate
    # end
    # fz_f = t-> charge_function(β, t, B, rate) - βf
    # method = Roots.A42()
    # lb = charged_β / rate - 1e-3
    # ub = charged_β/rate/0.2
    # t_opt = find_zero(fz_f, (lb, ub), method=method) 
    # return t_opt 
end

function charge_function_g(β::Real, t::Real, B::Real, rate::Real = g_default_rate)
    dβ = FiniteDiff.finite_difference_derivative(x->charge_function(x,t, B, rate), β)
    dt = FiniteDiff.finite_difference_derivative(x->charge_function(β,x, B, rate), t)
    return (dβ, dt)
end

function charge_function_h(β::Real, t::Real, B::Real, rate::Real = g_default_rate)
    f = (x)-> charge_function(x[1], x[2], B, rate)
    return FiniteDiff.finite_difference_hessian(f, [β, t])
end