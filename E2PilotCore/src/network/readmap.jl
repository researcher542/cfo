
# The global network dict to cache the network.
global g_net_dict = Dict{String, Network}()

function map_has_csv_data(region::String)
    if region == "oldmap" || occursin("test", region) || region == "eu" || occursin("eu.merged", region)
        return true
    end
    return false
end

"""
The region in a net is more like a name of map that identifies the map.

The continent is something like :us and :eu that represents a larger region
"""
function region2continent(region::String)
    if occursin("eu", region)
        return :eu
    else
        return :us
    end
end

"""
read from map data
bin: use the binary file or not; It seems 
use_saved: use the saved net in the global dict
"""
function readmapdata(region::String="test"; 
        isreduce::Bool = true, 
        bin::Bool=false, 
        use_saved::Bool=true,
        osm_flag::Bool=false,
        )::Network
    if haskey(g_net_dict, region) && use_saved
        return g_net_dict[region]
    end
    if bin
        (flag, net) = load_bin_map(region)
        if flag
            g_net_dict[region] = net
            return net
        end
    end
	if map_has_csv_data(region) && !osm_flag
        net = readoldmap(region)
    else
        net = Network(region)	
        @info "reading nodes"
        #readnodes_csv!(net, region)
        readnodes_osm!(net, region)

        @info "reading ways"
        readways_osm!(net, region)
	end
    net.continent = region2continent(region)

	n_edge = ne(net.graph)
    addways!(net)
    # check_self_loop(net)
    # @debug "after merging."
	if isreduce
		reduce_edge!(net)
	end
	add_junction!(net)
	n_edge_red = ne(net)

	@debug "Got $(nv(net)) nodes, $(n_junc(net)) junctions, $(nway(net)) ways, $(n_edge) edges, $(n_edge_red) edges after reduction."
	# @debug "$(n_junc(net)) junctions got."
	# @debug "$(nway(net)) ways got."
	# @debug "$(n_edge) total edges"
	# @debug "$(n_edge_red) edges after merge and reduction."
    g_net_dict[region] = net
    if bin save_bin_map(net, region) end
	return net
end


"""Load the binary map"""
function load_bin_map(region::String)
    @debug "loading binary map for $(region)..."
    net_data_path = joinpath(g_bin_net_data_dir, "$region.jld2") 
    if isfile(net_data_path)
        net = JLD2.load(net_data_path)["net"]
        @debug "Binary map for $(region) loaded"
        return (true, net)
    else
        @debug "Binary map for $(region) not exist."
        return (false, Network(""))
    end
end

"""Save the binary map"""
function save_bin_map(net::Network, region::String)
    @debug "saving binary map for $(region)..."
    net_data_path = joinpath(g_bin_net_data_dir, "$region.jld2") 
    JLD2.save(net_data_path, "net", net)
end