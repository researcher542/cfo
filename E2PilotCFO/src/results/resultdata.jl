

@with_kw mutable struct CfoResult
    alg::String = ""
    time::Float64 = NaN # total traveling time
    cost::Float64 = NaN # total cost
    src::Int = 0
    des::Int = 0
    idx::Int = 0 # the index in the faf data, should not be used, I guess...
    γ::Float64 = 0.0
    start_time::DateTime = DateTime(2022,1,1) 
    t_vec::Vector{CsTVar{Float64}} = []
    path::Vector{Int} = []
    β_vec::Vector{Float64} = []
    cap::Float64 = NaN # The battery capacity.
    β0::Float64 = NaN # The initial 
    meta::Dict{Symbol, Any} = Dict(
        :renewable_mul => 1.0,
    )
end

isvalid(a::Any) = false

function isvalid(res::CfoResult)
    if haserror(res) 
        return false
    elseif isnan(res.time)
        return false
    end
    return true
end

function haserror(res::CfoResult)
    if haskey(res.meta, :error) || haskey(res.meta, :err)
        return true
    end
    return false
end

function Base.haskey(res::CfoResult, key::Symbol)
    if hasfield(CfoResult, key)
        return true
    elseif haskey(res.meta, key)
        return true
    else
        return false
    end
end

function Base.getproperty(res::CfoResult, name::Symbol)
    if name in fieldnames(CfoResult)
        return Base.getfield(res, name)
    elseif haskey(res.meta, name)
        return res.meta[name]
    elseif name == :renewable_mul && !haskey(res.meta, :renewable_mul)
        return 1.0
    elseif name == :infeasible_flag && !haskey(res.meta, :infeasible_flag) && occursin("ice", res.alg)
        return false
    else
        throw(ErrorException("CfoResult has no property $name"))
    end
end

function Base.setproperty!(res::CfoResult, name::Symbol, x)
    if name in fieldnames(CfoResult)
        return Base.setfield!(res, name, x)
    else
        return res.meta[name] = x
    end
end

function Base.:(==)(res::CfoResult, res0::CfoResult)
    sym_list = [:idx, :alg, :src, :des, :γ, :start_time, :cap]
    if res.alg == "paso"
        sym_list = vcat(sym_list, [:t_ratio])
    end
    for sym in sym_list
        val1 = getproperty(res, sym)
        val2 = getproperty(res0, sym)
        if val1 != val2
            return false
        end
    end
    return true
end


function is_ice_truck(res::CfoResult)
    if occursin("ice", res.alg) && !occursin("practice", res.alg)
        return true
    else
        return false
    end
end


"""
Get carbon cost from the result. The unit is kg.

The carbon cost also accounts for the carbon footprint used at the source. (Note that in the new set of simulation, )
"""
function get_carbon_cost(res::CfoResult, prob)
    if is_ice_truck(res)
        # @show "ice" res.alg
        return dieselenergy2carbon(res.cost)
    else
        # @show "e-truck" res.alg
        cost = res.cost 
        # return cost / 3.6e6 / g_charge_eff
        ta = get_price_ta(prob, res.src, PredictPerfect())
        delta_beta =  max(0, res.β0 - res.β_vec[end])
        p = get_price(ta, res.start_time, 0.0)[1]
        # # @show p
        # cost -= delta_beta * p
        return cost / 3.6e6 / g_charge_eff
    end
end

function get_energy_cost(res::CfoResult)
    if is_ice_truck(res)
        return res.cost
    else
        e_cost = res.e_cost
        # e_cost += max(0, res.β0 - res.β_vec[end])
        return e_cost / g_charge_eff
    end
end
