"""
This is a self-implemented algorithm for time-expanded and battery-expanded graph.

This one is similar to the files for stage-expanded graph. Not sure if we can combine those two
"""

@enum TbNodeType begin
    InCharge = 1 # the node is one go into a charge station
    WaitCharge = 2 # the node is one wait at the charge station
    OutCharge = 3 # the node is one go out of the charge station
    RealRoad = 4 # the node of the real road
end

"""
A node in the time and battery expanded graph
"""
@with_kw struct TBNode
    orig_idx::Int ## the idx of the original graph
    ext_idx::Int ## the 
    itime::Int
    ibattery::Int
    beta::Float64 = 0.0 # the battery level of the node
    tau::Float64 = 0.0 # the time level of the node
    type::TbNodeType
end

@with_kw mutable struct TimeBatteryExpNetwork{G, DT} <: AbsNet
    orig_g::G 
    n_t_state::Int
    n_b_state::Int
    w_dict::DT = OrderedDict{Tuple{Int, Int}, Float64}() # The dict to store the weight of edges in the time and battery expanded graph.
    nodesdata::OrderedDict{Int, TBNode} # store the nodes with idx as a key

    out_neigh_dict::OrderedDict{Int, Vector{Int}} = OrderedDict{Int, Vector{Int}}() # The dict to store the buffer of outneighbors in the BT expanded network.

    # virtual_des::Int = -1 # The virtual destination node in the time and battery expanded graph.

    orig_des::Int = -1 # The original destination node, used to identify the virtual destination node.
    
    max_index::Int = -1 # The maximum index of the expanded graph, 

    # out_neigh_dict::OrderedDict{Int, Vector{Int}} = OrderedDict{Int, Vector{Int}}() # The dict to store the buffer of outneighbors in the stage ext network.
end

function Base.getproperty(exp_net::TimeBatteryExpNetwork, key::Symbol) 
    if key == :virtual_des
        return exp_net.max_index + 1
    else
       return getfield(exp_net, key)
    end
end

Base.copy(exp_net::TimeBatteryExpNetwork) = TimeBatteryExpNetwork(
    orig_g = exp_net.orig_g,
    n_t_state = exp_net.n_t_state,
    n_b_state = exp_net.n_b_state,
    max_index = exp_net.max_index,
    w_dict = deepcopy(exp_net.w_dict),
    nodesdata = deepcopy(exp_net.nodesdata),
    out_neigh_dict = deepcopy(exp_net.out_neigh_dict),
)

"""
For time and battery expanded graph.
"""
@with_kw mutable struct TbOption
    eps_t::Float64 = 60.0  # in seconds
    eps_b::Float64 = 10 * 3.6e6  # in J, default is 10 kWh
    mem_limit::Float64 = 16 * 1024 # in MB
    nv_limit::Int = 5e7 # the limit of the number of edges in the expanded graph
    ne_limit::Int = 2e8 # the limit of the number of edges in the expanded graph
    t_mul::Float64 = 1.0 / 60.0 # the multiplier for time
    b_mul::Float64 = 1.0 / 3.6e6 # the multiplier for battery
end


function orig2exp_idx(orig_idx::Int, it::Int, ib::Int, n_t_state::Int, n_b_state::Int, nt::TbNodeType)
    @assert it > 0 && it <= n_t_state+1
    @assert ib > 0 && ib <= n_b_state+1


    exp_idx = orig_idx

    exp_idx = exp_idx * (n_b_state+2) + ib

    exp_idx = exp_idx * (n_t_state+2) + it

    n_type = length(instances(TbNodeType))

    exp_idx = exp_idx * n_type + Int(nt)


    return exp_idx
end

function add_nodes!(exp_net::TimeBatteryExpNetwork, prob::AbsCfoProb, option::TbOption)

    @info "Adding nodes for the time and battery expanded graph..."
    (; net) = prob
    (; n_b_state, n_t_state, nodesdata) = exp_net
    (; eps_t, eps_b) = option
    
    for orig_node in net.nodesdata
        check_ve_limit(exp_net, option)
        orig_idx = orig_node.idx
        node_type_vec = [RealRoad]

        ## If the node is a charging station, we need to add two nodes for the in and out charge.
        if net.cs_flag_vec[orig_idx]
            node_type_vec = [InCharge, WaitCharge, OutCharge, RealRoad]
        end

        for ib in 1:n_b_state
            for it in 1:n_t_state+1
                for nt in node_type_vec
                    ext_idx = orig2exp_idx(orig_idx, it, ib, n_t_state, n_b_state, nt)
                    node = TBNode(
                        orig_idx = orig_idx, 
                        ext_idx = ext_idx, 
                        itime = it, 
                        ibattery = ib, 
                        beta = ib * eps_b,
                        tau = it * eps_t,
                        type = nt)
                    @assert !haskey(exp_net.nodesdata, ext_idx) (@exfiltrate; "The node $ext_idx is already in the expanded graph!")
                    exp_net.nodesdata[ext_idx] = node
                    exp_net.max_index = max(exp_net.max_index, ext_idx)
                end
            end 
        end
    end
    @info "Construced TBNet with $(length(exp_net.nodesdata)) nodes."
    return exp_net
end

function Graphs.add_edge!(exp_net::TimeBatteryExpNetwork, u::Int, v::Int, w::Float64)
    if !haskey(exp_net.nodesdata, u) 
        @warn "The node $u is not in the expanded graph!"
        return
    end
    if !haskey(exp_net.nodesdata, v) 
        @warn "The node $v is not in the expanded graph!"
        return
    end

    if !haskey(exp_net.w_dict, (u, v))
        exp_net.w_dict[(u, v)] = w
    end

    if !haskey(exp_net.out_neigh_dict, u)
        exp_net.out_neigh_dict[u] = Vector{Int}()
    end
    push!(exp_net.out_neigh_dict[u], v)
end

function check_ve_limit(exp_net::TimeBatteryExpNetwork, option::TbOption)
    if length(exp_net.nodesdata) > option.nv_limit
        @error "The number of nodes in the expanded graph is too large! $(length(exp_net.nodesdata)) > $(option.nv_limit)"
        throw(ErrorException("The number of nodes in the expanded graph is too large! $(length(exp_net.nodesdata)) > $(option.nv_limit)"))
    end
    if length(exp_net.w_dict) > option.ne_limit
        @error "The number of edges in the expanded graph is too large! $(length(exp_net.w_dict)) > $(option.ne_limit)"
        throw(ErrorException("The number of edges in the expanded graph is too large! $(length(exp_net.w_dict)) > $(option.ne_limit)"))
    end
    
end

function add_edges!(exp_net::TimeBatteryExpNetwork, prob::AbsCfoProb, option::TbOption)
    @info "Adding edges for the time and battery expanded graph..."
    (; src, des, T, ev, net) = prob
    (; n_b_state, n_t_state) = exp_net
    B = ev.cap

    (; eps_t, eps_b) = option

    t_begin = time()

    t_start = time()
    inode = 0
    for ex_node in values(exp_net.nodesdata)
        inode += 1
        if time() - t_begin > 10.0
            @debug "Adding edges for the time and battery expanded graph... $(time() - t_start) seconds, ne=$(length(exp_net.w_dict)), nv=$(length(exp_net.nodesdata)) inode=$inode"
            check_ve_limit(exp_net, option)
            t_begin = time()
        end
        orig_idx = ex_node.orig_idx
        if ex_node.type == InCharge
            for itau in ex_node.itime+1:n_t_state
                tw = (itau - ex_node.itime) * eps_t
                if tw < g_min_wait_time || tw > g_max_wait_time
                    continue
                end

                next_node_idx = orig2exp_idx(orig_idx, itau, ex_node.ibattery, n_t_state, n_b_state, WaitCharge)
                if haskey(exp_net.nodesdata, next_node_idx)
                    w = 0.0
                    add_edge!(exp_net, ex_node.ext_idx, next_node_idx, w)
                end
            end
        elseif ex_node.type == WaitCharge
            ## Add the edges for the in charge node, the edges must be outcharge
            for ibattery in ex_node.ibattery+1:n_b_state
                b0 = ex_node.beta
                bf = ibattery * eps_b

                tc = cf_inv(b0, bf, B)
                itc = round(Int, tc / eps_t, RoundUp) 
                tc_int = itc * eps_t
                itime = ex_node.itime + itc

                if itime > n_t_state
                    continue
                end

                next_node_idx = orig2exp_idx(orig_idx, itime, ibattery, n_t_state, n_b_state, OutCharge)
               
                w = cfo_objective(prob, orig_idx, b0, tc_int, ex_node.tau, prob.predict_mode)
                if haskey(exp_net.nodesdata, next_node_idx)
                    add_edge!(exp_net, ex_node.ext_idx, next_node_idx, w)
                end
            end
        elseif ex_node.type == OutCharge || ex_node.type == RealRoad
            #
            for orig_out_idx in Graphs.outneighbors(net, orig_idx)
                min_t, max_t = minmax_t(net, orig_idx, orig_out_idx, ev)
                imin_t = round(Int, min_t / eps_t, RoundUp)
                imax_t = round(Int, max_t / eps_t, RoundUp)
                node_type_vec = [RealRoad]
                if net.cs_flag_vec[orig_out_idx]
                    node_type_vec = [InCharge, RealRoad]
                end
                for it in imin_t:imax_t
                    t = it * eps_t
                    e_cost = energy_cost_on_road(net, orig_idx, orig_out_idx, ev, t)
                    next_b = ex_node.beta - e_cost
                    ib = round(Int, next_b / eps_b, RoundDown)
                    itau = ex_node.itime + it
                    if itau > n_t_state || ib < 1 || ib > n_b_state
                        continue
                    end
                    for nt in node_type_vec
                        next_node_idx = orig2exp_idx(orig_out_idx, itau, ib, n_t_state, n_b_state, nt)
                        if haskey(exp_net.nodesdata, next_node_idx)
                            w = e_cost * 1e-9
                            add_edge!(exp_net, ex_node.ext_idx, next_node_idx, w)
                        end
                    end
                end
            end
        else
            @error "The node type is not supported!"
            return exp_net
        end
    end

    return exp_net
end

function check_memory_limit(option)
    used_mem = ep.get_cur_program_mem()
    @show used_mem
    if used_mem > option.mem_limit
        # Threads.terminate(task)
        @error "The memory limit is exceeded!"
        throw(ErrorException("The memory limit is exceeded! mem=$used_mem, limit=$(option.mem_limit)"))
    end
    return 
end

function construct_bte_net(prob::AbsCfoProb, option)
    (; src, des, T, ev, net) = prob
    B = ev.cap
    (; eps_t, eps_b) = option
    n_t_state = round(Int, T / eps_t, RoundUp) + 1
    n_b_state = round(Int, B / eps_b, RoundUp) + 1
    @info "Constructing TB Net with n_t_state = $n_t_state, n_b_state = $n_b_state "

    exp_net = TimeBatteryExpNetwork(
        orig_g = prob.net.g, 
        nodesdata = OrderedDict{Int, TBNode}(),
        w_dict = OrderedDict{Tuple{Int, Int}, Float64}(),
        out_neigh_dict = OrderedDict{Int, Vector{Int}}(),
        n_t_state = n_t_state, 
        n_b_state = n_b_state,
        )
    
    add_nodes!(exp_net, prob, option)
    check_ve_limit(exp_net, option)
    # check_memory_limit(option)
    add_edges!(exp_net, prob, option)
    # check_memory_limit(option)
    

    return exp_net
end

function Graphs.nv(exp_net::TimeBatteryExpNetwork)
    return exp_net.virtual_des + 1
end

function num_node(exp_net::TimeBatteryExpNetwork)
    return length(exp_net.nodesdata)
end

Graphs.ne(exp_net::TimeBatteryExpNetwork) = length(exp_net.w_dict)
function Graphs.outneighbors(exp_net::TimeBatteryExpNetwork, ex_idx::Int) 

    if !haskey(exp_net.out_neigh_dict, ex_idx)
        return []
    end

    orig_idx = exp_net.nodesdata[ex_idx].orig_idx

    if orig_idx == exp_net.orig_des
        return exp_net.virtual_des
    end

    return exp_net.out_neigh_dict[ex_idx]
end

function ex_get_dist(net::TimeBatteryExpNetwork, u::Int, v::Int, lam, prob)
    if v == net.virtual_des
        ## We need to account for the initial battery level here
        init_ci = get_init_battery_ci(prob, prob.predict_mode)
        b0 = prob.β0
        cost =  init_ci * (b0 - net.nodesdata[u].beta)
        return cost
    end
    if haskey(net.w_dict, (u, v))
        return net.w_dict[(u, v)]
    end
    return Inf 
end

 function node_vec2sol(prob, nd_vec::Vector{TBNode})
    @info "node_vec2sol" nd_vec
    tc = 0.0 
    tw = 0.0
    beta0 = prob.β0
    if nd_vec[1].type == InCharge
        path = [nd.orig_idx for nd in nd_vec[3:end]]
        t_vec = [nd_vec[i+1].tau - nd_vec[i].tau for i in 3:length(nd_vec)-1]
        tw = nd_vec[2].tau - nd_vec[1].tau
        tc = nd_vec[3].tau - nd_vec[2].tau
        beta0 = nd_vec[1].beta
    else
        ## I guess the this is for the first path where we do not charge at the beginning.
        path = [nd.orig_idx for nd in nd_vec[1:end]]
        t_vec = [nd_vec[i+1].tau - nd_vec[i].tau for i in 1:length(nd_vec)-1]
    end
    sol = DsSubsolution(
        path = path,
        t_vec = t_vec,
        tc = tc,
        tw = tw,
        τ = nd_vec[1].tau, 
        β = beta0,
    ) 
    return sol
end

"""
Run the shortest path algorithm on the time and battery expanded graph. 
"""
function cfo_time_battery_expanded(prob0::AbsCfoProb, option, exp_net0 = nothing)
    prob = get_cfo_prob_fix(prob0) 
    if isnothing(exp_net0) 
        exp_net = construct_bte_net(prob, option)
    else
        exp_net = exp_net0
    end
    exp_net.orig_des = prob.des

    ib = round(Int, prob.β0 / option.eps_b, RoundDown)
    ex_src = orig2exp_idx(prob.src, 1, ib, exp_net.n_t_state, exp_net.n_b_state, RealRoad)

    task = Threads.@spawn begin
        path = ep.shortest_path(Astar(), exp_net, ex_src, exp_net.virtual_des
            ;
            getdist = ex_get_dist,
            breakearly = false,
            visitonce = false,
            getdistargs = (prob, ),
        )
        return path
    end

    errormonitor(task)
    used_mem = ep.get_cur_program_mem()
    @show used_mem
    while !Threads.istaskdone(task)
        sleep(1.0)
        # check_memory_limit(option)
    end

    path = Threads.fetch(task)
    

    node_vec = [
        exp_net.nodesdata[path[i]] for i in 1:length(path)-1
    ]
    cost = sum(
        ex_get_dist(exp_net, path[i], path[i+1], 0.0, prob) for i in 1:length(path)-1)
    @show cost / 3.6e6

    sol_vec = DsSubsolution[]
    inner_node_vec = TBNode[]
    
    for node in node_vec
        if node.type == InCharge
            # if the node is in charge, we need to start a new subsolution
            sol = node_vec2sol(prob, vcat(inner_node_vec, node))
            push!(sol_vec, sol)
            inner_node_vec = [node]
        else
            push!(inner_node_vec, node)
        end
    end
    if length(inner_node_vec) > 1
        sol = node_vec2sol(prob, inner_node_vec)
        push!(sol_vec, sol)
    end
    # convert the node vector to the primal solution.
    primal = DsPrimal(; sub_sol_vec = sol_vec)

    return primal
end

"""
Check if the forward and backward index mapping is correct.
"""
function check_idxmap_exp_net(exp_net::TimeBatteryExpNetwork)
    @info "checking the index mapping of the expanded network..."
    for (i, nd) in exp_net.nodesdata
        orig_idx = nd.orig_idx
        it = nd.itime
        ib = nd.ibattery
        nt = nd.type
        ex_idx = orig2exp_idx(orig_idx, it, ib, exp_net.n_t_state, exp_net.n_b_state, nt)
        if ex_idx != i
            @error "The index mapping is not correct! $i != $ex_idx"
            return false
        end
    end
    return true
end

function report_failed_g_score_on_sp(exp_net, g_score, orig_sp)
    (;n_b_state, n_t_state) = exp_net
    for orig_idx in orig_sp
        @info "orig_idx: $(orig_idx)"
        break_flag = false
        for it in 1:n_t_state
            for ib in n_b_state:-1:1
                # @info "ib: $(ib), it: $(it)"
                for nt in instances(TbNodeType)
                    ext_idx = orig2exp_idx(orig_idx, it, ib, n_t_state, n_b_state, nt)
                    if !isinf(g_score[ext_idx])
                        @info "orig_idx: $(orig_idx), ib: $(ib), it: $(it), nt: $(nt), ext_idx: $(ext_idx), g_score: $(g_score[ext_idx])"
                        break_flag = true
                        break
                    end
                end
                break_flag && break
            end
            break_flag && break
        end
    end
    
end