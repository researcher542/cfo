"""
Read the old data map from Wenjie's code.
@Note: Need to be careful about the indexing, since julia's index starts from 1.
"""

"""
Read the old data map from Wenjie's code.
or read data from the test folder
"""
function readoldmap(name::String = "oldmap")::Network
	net = Network(name)
	@info "reading old map data: $name"
	datadir = joinpath(k_root_dir, "data", "map", name)

    @info "reading old map node"
	readoldmap_node!(net, datadir)

    @info "reading old map edge"
	readoldmap_edge!(net, datadir)

	return net
end

"""
"""
function export_oldmap_node(net::Network, datadir::String)
	datapath = joinpath(datadir, "node_new.csv")
    @info "exporting node data to $datapath"

    node_vec = net.nodesdata
    df = DataFrame(
        name = String[nd.name for nd in node_vec],
        lat = Float64[nd.lat for nd in node_vec],
        lon = Float64[nd.lon for nd in node_vec],
        ele = Float64[nd.ele for nd in node_vec],
        id = Int64[nd.id for nd in node_vec],
    ) 
    CSV.write(datapath, df)

end

function readoldmap_node!(net::Network, datadir::String)
	datapath = joinpath(datadir, "node.csv")
	df = DataFrame(CSV.File(datapath))
    cnt = 1
    # n_node = length(df.id)
	for row in eachrow(df)
		id = row.id
        name = ismissing(row.name) ? Symbol() : Symbol(row.name)
		# node = Node(row.lat, row.lon, row.ele, id, cnt, name)
		node = StaticNode(row.lat, row.lon, row.ele, id, cnt, name)
		addnode!(net, node)
        cnt += 1
        if cnt % 1000_000 == 0
            @debug "processing $(cnt)-th node"
        end
	end
	return net
end

function isoldmap_region(region_sym::Symbol)
    return  region_sym == :oldmap || region_sym == :testsmall || region_sym == :test4node
end

function readoldmap_edge!(net::Network, datadir::String)
	datapath = joinpath(datadir, "edge.csv")
	df = DataFrame(CSV.File(datapath))
    cnt = 1
    region_sym = Symbol(net.region)
    src::Int64 = -1
    des::Int64 = -1
    speeddata_vec = Pair{Tuple{Int, Int}, Float64}[]
	for row in eachrow(df)
        src = row.src
        des = row.des
        (src == des) && continue
        # @show src des
        if isoldmap_region(region_sym)
		    src = src + 1
		    des = des + 1
        end
		add_edge!(net.graph, src, des)
        # net.speeddata[(src,des)] = mph2ms(row.speed_pos)
        # set_speeddata!(net, src, des, mph2ms(row.speed_pos))
        push!(speeddata_vec, (src, des) => mph2ms(row.speed_pos))
        if isoldmap_region(region_sym)
		    add_edge!(net.graph, des, src)
            push!(speeddata_vec, (des, src) => mph2ms(row.speed_neg))
            # set_speeddata!(net, des, src, mph2ms(row.speed_neg))
            # net.speeddata[(des,src)] = mph2ms(row.speed_neg)
        end
        if cnt % 1000_000 == 0
            @debug "processing $(cnt)-th edge"
        end
        cnt += 1
        
	end
    # @show length(speeddata_vec)
    @time d = Dict(speeddata_vec)
    net.speeddata = d
    return net
end

set_speeddata!(net, src::Int, des::Int, speed::Real) = net.speeddata[(src, des)] = speed