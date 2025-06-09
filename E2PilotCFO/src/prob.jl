

abstract type AbsCfoProb end

"""
This struct stores the large data to deine a CfoProb. 
We should avoid copying the data in this struct to save memory.
"""
@with_kw struct CfoProbData{S<:EV, CT <: AbstractDict, CSDT <: AbstractDict}
    net::Network
    ev::S
    carbon_dataset::AbstractCarbonDataset = CambiumDataset()
    carbon_dict::CT = get_carbon_dict(carbon_dataset) # the carbon dict that stores the information of carbon-intensity.
    carbon_predict_dict::CT = get_carbon_dict_predict(carbon_dataset) # the carbon dict that stores the information of carbon-intensity.
    cs_dict::CSDT
end

@with_kw mutable struct CfoProb{DT<:Ref} <: AbsCfoProb
    fix_data::DT
    src::Int
    des::Int
    β0::Float64
    T::Float64
    γ::Float64 = 0.0
    start_time::DateTime
    objtype::AbstractObjective = ObjCarbon()
    scenario::String = g_default_scenario
    predict_mode::AbstractPredictionMode = PredictPerfect()
    odset::AbstractODSet
end

"""
A fixed type of CfoProb that stores the data in the struct. 
We use this for type stability.
"""
@with_kw struct CfoProbFix{DT<:Ref, OT <: AbstractObjective, PMT <: AbstractPredictionMode} <: AbsCfoProb
    fix_data::DT
    src::Int
    des::Int
    β0::Float64
    T::Float64
    γ::Float64 = 0.0
    start_time::DateTime
    objtype::OT = ObjCarbon()
    scenario::String = g_default_scenario
    predict_mode::PMT = PredictPerfect()
end

function get_cfo_prob_fix(p::AbsCfoProb)
    return CfoProbFix(
        fix_data = p.fix_data,
        src = p.src,
        des = p.des,
        β0 = p.β0,
        T = p.T,
        γ = p.γ,
        start_time = p.start_time,
        objtype = p.objtype,
        scenario = p.scenario,
        predict_mode = p.predict_mode
    )
end

function Base.getproperty(p::AbsCfoProb, sym::Symbol)
    if hasfield(CfoProb, sym)
        return getfield(p, sym)
    else
        return getproperty(p.fix_data[], sym)
    end
    
end

function Base.copy(p::AbsCfoProb)
    return CfoProb(;
        fix_data = p.fix_data,
        src = p.src,
        des = p.des,
        β0 = p.β0,
        T = p.T,
        γ = p.γ,
        start_time = p.start_time,
        objtype = p.objtype,
        scenario = p.scenario,
        predict_mode = p.predict_mode,
        odset = p.odset,
    )
end


# end