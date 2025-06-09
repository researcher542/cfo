
abstract type AbsDictWrapper end

struct DictWapper{T} <: AbsDictWrapper where T <: AbstractDict
    d::T
end

function Base.getproperty(dw::AbsDictWrapper, sym::Symbol)
    if sym in fieldnames(typeof(dw))
        return getfield(dw, sym)
    else
        res = get(dw.d, sym, nothing)
        if isnothing(res)
            error("Dict has not key $sym")
        end
        return res
    end
end