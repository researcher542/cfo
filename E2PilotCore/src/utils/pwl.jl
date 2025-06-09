"""
Utilities for Piecewise linear functions.
"""

"""
"""
struct PWL
    x::Vector{Float64}
    y::Vector{Float64}
end

function PWL(; x::Vector{Float64} = Float64[], y::Vector{Float64} = Float64[])
    @assert length(x) == length(y)
    @assert length(x)>=2
    if !issorted(x)
        perm::Vector{Int} = sortperm(x)
        return PWL(x[perm], y[perm])
    end
    @assert issorted(x)
    return PWL(x,y) 
end

"""
Given y, return the x value
"""
function inv(pwl::PWL, y::Real; ext::Bool = false)
    pwl_inv = PWL(;x=pwl.y, y=pwl.x) 
    return pwl_inv(y, ext)
end

function interpolate_pwl(x::T, x_vec::Vector{T}, y_vec::Vector{T}, ext::Bool = false) where T <: Real
    # If it is not in the domain.
    tol = 1e-10
    if !ext && (x < (first(x_vec) - tol) || x > (last(x_vec) + tol) )
        @warn "get a value out of domain"
        error("Get a value out of domain")
        return NaN
    end
    if x < first(x_vec)
        x1,x2 = x_vec[1], x_vec[2]
        y1,y2 = y_vec[1], y_vec[2]
    elseif x > last(x_vec)
        x1,x2 = x_vec[end-1], x_vec[end]
        y1,y2 = y_vec[end-1], y_vec[end]
    elseif x == first(x_vec)
        return first(y_vec)
    elseif x == last(x_vec)
        return last(y_vec)
    else
        idx = searchsortedfirst(x_vec, x)
        if idx == length(x_vec) + 1
            @show x x_vec
        end
        x1,x2 = x_vec[idx-1], x_vec[idx]
        y1,y2 = y_vec[idx-1], y_vec[idx]
    end
    (x1 == x) && return y1
    (x2 == x) && return y2
    deri = (x2==x1) ? 0.0 : (y2-y1)/(x2-x1)
    return deri*(x-x1)+y1
    
end

"""
ext: if it is ok to evualate the value out of x_value
"""
function (pwl::PWL)(x::Real, ext::Bool = false)
    interpolate_pwl(x, pwl.x, pwl.y, ext) 
end

minmax_x(pwl::PWL) = (first(pwl.x), last(pwl.x))

