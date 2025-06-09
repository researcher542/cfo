"""
Copy the corresponding file from Optim.jl for custimized function design
"""

"""
## Description
The `GoldenSection` method seeks to minimize a univariate function on an interval
`[a, b]`. At all times the algorithm maintains a tuple of three minimizer candidates
`(c, d, e)` where ``c<d<e`` such that the ratio of the largest to the smallest interval is the Golden Ratio.

## References
https://en.wikipedia.org/wiki/Golden-section_search

return: (minimum, minimizer)
"""
function golden_section(f::Function, x_lower::T, x_upper::T;
     rel_tol::T = sqrt(eps(T)),
     abs_tol::T = eps(T),
     iterations::Int = 1_000,
     show_trace::Bool = false,
    ) where T <: AbstractFloat

    if x_lower > x_upper
        error("x_lower must be less than x_upper")
    end
    golden_ratio::T = 0.5 * (3.0 - sqrt(5.0)) # 1/phi^2 where phi=1.618

    initial_lower = x_lower
    initial_upper = x_upper

    new_minimizer = x_lower + golden_ratio*(x_upper-x_lower)
    new_minimum = f(new_minimizer)
    best_bound = "initial"
    f_calls = 1 # Number of calls to f

    iteration = 0
    converged = false

    while iteration < iterations

        x_tol = rel_tol * abs(new_minimizer) + abs_tol

        x_midpoint = (x_upper+x_lower)/2

        if abs(new_minimizer - x_midpoint) <= 2*x_tol - (x_upper-x_lower)/2
            converged = true
            break
        end

        iteration += 1

        if x_upper - new_minimizer > new_minimizer - x_lower
            new_x = new_minimizer + golden_ratio*(x_upper - new_minimizer)
            new_f = f(new_x)
            f_calls += 1
            if new_f < new_minimum
                x_lower = new_minimizer
                best_bound = "lower"
                new_minimizer = new_x
                new_minimum = new_f
            else
                x_upper = new_x
                best_bound = "upper"
            end
        else
            new_x = new_minimizer - golden_ratio*(new_minimizer - x_lower)
            new_f = f(new_x)
            f_calls += 1
            if new_f < new_minimum
                x_upper = new_minimizer
                best_bound = "upper"
                new_minimizer = new_x
                new_minimum = new_f
            else
                x_lower = new_x
                best_bound = "lower"
            end
        end

        if show_trace
            @debug "iteration $iteration" x_lower x_upper new_minimizer new_minimum
        end
    end

    return (new_minimum,new_minimizer)
end