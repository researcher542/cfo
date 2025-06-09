
str2int(s::String) =parse(Int64, s)
str2float(s::String) = parse(Float64, s)

function Base.hash(tp::Tuple{Int,Int}, h::UInt)
    hash(tp[1], hash(tp[2], h)) 
end


pos(v::Vector{T}) where T <: Real = pos.(v)
"""
Return max(0, x)
"""
pos(x::Real) = (x > 0.0 ? x : 0.0)

"""
This is tricky, the built-in can be slow. not sure why...
"""
function Base.min(x1::Float64, x2::Float64)
    return ifelse(x1<=x2, x1, x2)
    # return x1 <= x2 ? x1 : x2
end

remove_nan(vec::Vector{Float64}) = vec[findall(!isnan, vec)]

"""
The enum has the naming convention kENUMNAME.
This function will convert it to lower-case without the prefix "k"
"""
function enum2str(e)
	str = "$(e)"
	str = lowercase(str[2:end])
	return str
end

function str2enum(str::String)
	str1 = "k"*uppercase(str)
	sym = Symbol(str1)
	return eval(sym)
end

function searchsortedlast(v, x::Float64)
 	# theta_idx = searchsortedlast(theta_range, theta, 1, n, Base.Forward)
	# theta_idx = searchsortedlast(theta_range, theta)   
    n::Int = length(v)
    lo::Int = 1
    hi::Int = n + 1
    @inbounds while lo < hi - 1
        m::Int = Base.Sort.midpoint(lo, hi)
        val::Float64 = v[m]
        if x < val
            hi = m
        else
            lo = m
        end
    end
    return lo
end



"""convert theta to index"""
function index_in_range(range::StepRangeLen, val::Real)
    step = range[2] - range[1]
    idx1 = (val - range[1])/step
    # THis is tricky....
    idx = Integer(round(idx1, RoundUp)) + 1
    #@show step, idx1, idx, val
    if idx == 0
        idx = 1
    elseif idx > length(range)
        idx = length(range)
    end
    return idx
end

"""
The self-defined find_zero, to handle the boundary case. 
The function should be increasing.
return: the root.

--- 
keyword arguments: xatol, xrtol, atol, rtol
"""
function find_zero(f::Function, bnds::Tuple, method = Roots.A42(); kwargs...)
    lb = bnds[1]
    ub = bnds[2]
    fval_min = f(lb)
    fval_max = f(ub)
    if fval_min >= 0.0
        return lb
    elseif fval_max <= 0.0
        return ub
    end
    @assert fval_min <= fval_max "fval_min=$fval_min fval_max=$fval_max"
    # @show fval_max fval_min bnds kwargs
    return Roots.find_zero(f, bnds, method; kwargs...)
end


function findnearest(A::Vector{T}, t::T) where T 
    idx::Int = -1
    min_dis = Inf
    near_val = Inf
    @inbounds for (i, val) in enumerate(A)
        dis = abs(val-t)
        if min_dis > dis
            idx = i
            near_val = val 
        end
    end
    return idx,near_val
    # findmin(abs.(A.-t))[2]
end

"""
copy it from Base.task.jl, may be changed in the future?
"""
function set_task_tid(t::Task, tid::Int)
    # This may changed in 1.9, should be careful with this.
    @assert(VERSION <= v"1.10")
    @assert(1<= tid <= Threads.nthreads())
    ccall(:jl_set_task_tid, Cint, (Any, Cint), t, tid-1)
end

function getrange(start::Float64, stop::Float64, step::Float64)
    l::Int = Integer(round((stop-start)/step, RoundUp)) + 1
    if stop == start
        l = 1 
    end
    # r = range(start, stop, l)
    r = LinRange(start, stop, l)
    # vec::Vector{Float64} = collect(range(start, stop, l))
    return r
end

Base.copy(x::T) where T = T([getfield(x, k) for k in fieldnames(T)]...)

"""
"""
function is_in_hpc() 
    if occursin("junyansu2", homedir())
        if !is_sdsc_server()
            return true
        end
    end
    return false
end

is_cuhk_ie_server() = false

function is_sdsc_server() 
    return occursin("junyansu2", homedir()) && isfile(joinpath(homedir(), "issdsc"))
end

is_in_desktop() = occursin("sujy", homedir())

"""
Get the memory used by the current program.

Return the memory usage in MB.
"""
function get_cur_program_mem()
    pid = getpid()   
    output = read(`ps -p $pid -o %mem,rss,vsz`, String)
    # Parse the output
    lines = split(output, '\n')  # Split output into lines
    header = lines[1]           # First line (header)
    data = lines[2]             # Second line (actual data)

    # Extract memory info
    fields = split(strip(data))  # Split the data line into fields
    percent_mem = fields[1]      # First field: %MEM
    rss_kb = parse(Int, fields[2])  # Second field: RSS in kilobytes
    vsz_kb = parse(Int, fields[3])  # Third field: VSZ in kilobytes

    # Convert RSS to MB
    rss_mb = rss_kb / 1024
    vsz_mb = vsz_kb / 1024

    # Print the results
    # println("Memory Usage for Julia Process:")
    # println("  %MEM: $percent_mem%")
    # println("  RSS: $rss_mb MB")
    # return vsz_mb
    return rss_mb
end
