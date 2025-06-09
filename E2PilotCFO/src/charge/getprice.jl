

function get_price_ta(prob::AbsCfoProb, idx::Int, ci_predictmode::AbstractPredictionMode) 
    carbon_dict = prob.carbon_dict
    if ci_predictmode != PredictPerfect()
        carbon_dict = prob.carbon_predict_dict
    end
    return _get_price_ta(prob.fix_data[], carbon_dict, idx, prob.scenario)
end

"""
A function barrier
"""
_get_price_ta(fix_data::CfoProbData, carbon_dict, idx::Int, scenario::String) = _get_price_ta(fix_data.cs_dict, carbon_dict, idx, scenario)

function _get_price_ta(cs_dict, carbon_dict, i::Int, scenario::String)
    cs::ChargeStation = cs_dict[i]
    ta = carbon_dict[scenario][cs.region]    
    return ta
end

function _get_price_ta(net::Network, carbon_dict, i::Int, scenario::String)
    nd = getnode(net, i)
    state = latlon2subregion(nd.lat, nd.lon, :us)
    reg = g_us_state_abbr_dict[state]
    ta = carbon_dict[scenario][reg]    
    return ta
end


"""
st: the start_time 
return the price and the derivative of price w.r.t. τ, in price/s
"""
function get_price(price_ta::TS.TimeArray, st::DateTime, τ::Real)
    @assert τ >= 0 
    d_t = st + τ2second(τ)
    # ari_t rounded to hour
    ari_t_h = round(d_t, Hour, RoundDown)

    ts1 = TS.timestamp(price_ta)[1]
    if ts1 > st
        @warn "start time is earlier than the price array" ts1 st
        @assert false
    end


    # we use a quadratic function to approximate the function, 
    # this is to make the price function smoother.
    # d_ means the type is DateTime
    # d_t0 = ari_t_h - Hour(1)
    d_t1 = ari_t_h
    d_t2 = ari_t_h + Hour(1)
    # price0 = unsafe_getindex(price_ta, d_t0)
    price1::Float64 = unsafe_getindex(price_ta, d_t1)
    price2::Float64 = unsafe_getindex(price_ta, d_t2)
    # convert DateTime to milli-seconds
    # t = (d_t-st).value
    # t0, t1, t2 = (d_t0-st).value, (d_t1-st).value, (d_t2-st).value
    # Lagrange basis function
    # price = (t-t1)*(t-t2)/((t0-t1)*(t0-t2)) * price0 
    # price += (t-t0)*(t-t2)/((t1-t0)*(t1-t2)) * price1
    # price += (t-t0)*(t-t1)/((t2-t0)*(t2-t1)) * price2

    price_d = Inf

    # The difference between two DateTime is in ms, change it to hours
    ratio = (st + Second(round(τ)) - ari_t_h).value / 1000.0 / 3600.0
    # @debug ratio price1 price2
    price = (1-ratio)*price1 + ratio*price2
    # price_d = (price2-price1)/3600.0
    return price::Float64, price_d::Float64
end




# function unsafe_getindex_fast(ta::TS.TimeArray, time::DateTime, debug::Bool = false)
#     ts = TS.timestamp(ta)
# 
#     ts1 = ts[1]
#     ts2 = ts[2]
# 
#     val_v = TS.values(ta)
# 
#     idx = 1 + round(Int, (time - ts1).value / (ts2 - ts1).value) 
# 
#     if idx == length(ts) + 1
#         # @warn "time not found." time maxlog=10
#         @warn "time not found." time 
#         idx -= 1
#     end
# 
#     return val_v[idx]
# end

"""
The getindex function in TimeSeries package is very slow since it create a lot of array...

ignore_day: if true, ignore the day, only consider the time. Note that in NERL data, there is only month-hour data. So we need to ignore the day information.
"""
function unsafe_getindex(ta::TS.TimeArray, time0::DateTime, ignore_day::Bool = false, debug::Bool = false)
    time = time0
    if ignore_day
        time = DateTime(Year(time).value, Month(time).value, 1, Hour(time).value, Minute(time).value, Second(time).value)
    end
    ts = TS.timestamp(ta)
    val_v = TS.values(ta)
    idx = searchsortedfirst(ts, time)
    if idx == length(ts) + 1
        # @warn "time not found." time maxlog=10
        @warn "time not found." time 
        idx -= 1
    end
    # debug && @show idx
    return val_v[idx]
end

τ2second(τ) = Second(round(τ))


