

"""
The state to store the previously called values, such that we can reuse them and speed up.
"""
@with_kw mutable struct CfoObjState{TA}
    cs_idx::Int = -1
    beta::Float64 = -1.0
    tc::Float64 = -1.0
    tau::Float64 = -1.0
    price::Float64 = -1.0
    obj::Float64 = -1.0
    price_ta::TA = TS.TimeArray([DateTime(1)], [1.0])
end


# function cfo_objective(prob::AbsCfoProb, idx::Int, β0::Float64, tc::Float64, τ::Float64)
#     obj::Float64 = cfo_objective_low(prob, idx, β0, tc, τ, prob.objtype)
#     return obj
# end


"""
i: the index in the road segment network
"""
function cfo_objective(prob::AbsCfoProb, idx::Int, β0::Float64, tc::Float64, τ::Float64, ci_predictmode::AbstractPredictionMode, s::Union{CfoObjState, Missing} = missing)
    if prob.objtype == ObjCarbon()
        obj::Float64 = cfo_objective_low(prob, idx, β0, tc, τ, prob.objtype, ci_predictmode, s)
    else
        obj = cfo_objective_low(prob, idx, β0, tc, τ, prob.objtype)
    end
    return obj
end


"""
cfo objective with multiple dispatch

cfo_obj_state: the state to store the previously called values, such that we can reuse them and speed up.
"""
function cfo_objective_low(prob::AbsCfoProb, idx::Int, β::Real, t::Real, τ::Real, objtype::ObjCarbon, ci_predictmode::AbstractPredictionMode, cfo_obj_state::Union{Missing, CfoObjState})
    (;start_time::DateTime, ) = prob
    # ev = prob.fix_data[].ev
    # fix_data = get_fix_data(prob)
    B = prob.fix_data[].ev.cap
    Δβ::Float64 = (charge_function(β, t, B) - β)

    # price_ta = get_price_ta(prob, idx, ci_predictmode)

    # price::Float64, _ = get_price(price_ta, start_time, τ) 

    ###########
    if !ismissing(cfo_obj_state) && (idx == cfo_obj_state.cs_idx) 
        price_ta = cfo_obj_state.price_ta
    else
        price_ta = get_price_ta(prob, idx, ci_predictmode)
        if !ismissing(cfo_obj_state)
            cfo_obj_state.price_ta = price_ta
        end
    end

    if !ismissing(cfo_obj_state) && (idx == cfo_obj_state.cs_idx) && (τ == cfo_obj_state.tau)
        price = cfo_obj_state.price
    else
        price::Float64, _ = get_price(price_ta, start_time, τ) 
    end
    ###########3

    # price::Float64, _ = get_price(price_ta, start_time, τ) 
    # @debug "used price" price τ idx get_price(price_ta, start_time, 0.0)
    obj::Float64 = price*Δβ

    ## update the state
    if !ismissing(cfo_obj_state)
        cfo_obj_state.cs_idx = idx
        cfo_obj_state.beta = β
        cfo_obj_state.tc = t
        cfo_obj_state.tau = τ
        cfo_obj_state.price = price
        cfo_obj_state.obj = obj
    end
    #########

    return obj
end

function cfo_objective_low(prob::AbsCfoProb, idx::Int, β, t, τ, type::ObjEnergy)
    B = prob.ev.cap
    Δβ::Float64 = (charge_function(β, t, B) - β)
    return Δβ
end

function cfo_objective_low(prob::AbsCfoProb, idx::Int, β, t, τ, type::ObjTime)
    return t 
end




