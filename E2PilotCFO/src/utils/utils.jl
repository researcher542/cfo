
include("sim.jl")
include("shortestpath_all.jl")
include("oddata.jl")


function cs_vec2cs_flag_vec(net::Network, cs_vec::Vector{ChargeStation})
    n_node = nv(net)
    cs_idx_vec = [cs.idx for cs in cs_vec]
    cs_flag_vec = zeros(Bool, n_node)
    cs_flag_vec[cs_idx_vec] .= true
    return cs_flag_vec
end


function get_next_b(net::Network, i::Int, j::Int, t::Real, β::Real, ev::EV; no_clamp::Bool =false)
    Δb = Δsoc(net, i, j, t, β, ev)  
    next_b = β - Δb
    if no_clamp
        return next_b
    else
        return min(ev.cap, next_b)
    end
    #return next_b
end

function print_sparse(val::Float64, prefix)
    @debug prefix,val
end

function print_sparse(mat::SparseMatrixCSC, prefix)
    for (x,y,v) in zip(findnz(mat)...)
        @debug prefix,(x,y), v
    end
end

function print_sparse(vec::Vector, prefix)
    for (i,val) in enumerate(vec)
        if val != 0.0
            @debug prefix, i, val
        end
    end
end

function relative_error(x1, x2)
    if x1 == 0 && x2 == 0
        return 0.0
    end
    if x2 == 0 || x1 == 0
        return abs(x1-x2)
    end
    max_x = max(abs(x1), abs(x2))
    rel_err0 =  abs((x1-x2)/max_x)
    abs_err = max_x
    return min(abs_err, rel_err0)
        # return rel_err0
end



