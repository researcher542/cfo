"""
"""


@with_kw mutable struct ChargeStation
    lat::Float64 = Inf
    lon::Float64 = Inf
    idx::Int64 = -1 # the index of node in the network
    region::String = ""
    min_price::Float64 = 0.0
    max_price::Float64 = 0.0
    avg_price::Float64 = 0.0
    idx_r::Int64 = -1 # the index of node of corresponding real node
    # is_virtual::Bool = false
end

"""
Add the node of charge station to the network without separate the node.
"""
function add_charge_station_no_sep!(net::Network, cs_vec::Vector{ChargeStation})
    net.cs_vec = cs_vec
    net.cs_flag_vec = cs_vec2cs_flag_vec(net, cs_vec)
    net.cs_dict = OrderedDict([cs.idx for cs in cs_vec] .=> cs_vec)
    # add_junction!(net; cs_net_flag=true)
    return net
end

"""
Add the node of charge station to the network
For a path 1->2->3, when need to add two nodes to model one charge station at (2): 1->2->3 && 1->2->4->5->3, where 2,4,5 has the same location.
2->4 is the waiting edge, 4->5 is the charing edge, and 5->3 is a real road
"""
function add_charge_station!(net::Network, cs_vec::Vector{ChargeStation})
    cs_idx = nv(net) 
    net.cs_start_idx = cs_idx + 1
    if !net.sep_cs_node_flag
        return add_charge_station_no_sep!(net, cs_vec)
    end


    for cs in cs_vec 
        wait_idx = cs_idx + 1
        real_nd_idx = cs.idx
        node1 = getnode(net, real_nd_idx)
        wait_node = deepcopy(node1)
        wait_node.idx = wait_idx

        # tricky...
        cs.idx = wait_idx
        cs_node = deepcopy(node1)
        cs_idx += 2
        cs_node.idx = cs_idx

        addnode!(net, wait_node)
        addnode!(net, cs_node)
        for nd in outneighbors(net, real_nd_idx)
            spd = get_max_speed(net, real_nd_idx, nd)
            addway!(net, cs_idx, nd)
            net.speeddata[(cs_idx, nd)] = spd
        end
        addway!(net, real_nd_idx, wait_idx)
        addway!(net, wait_idx, cs_idx)
    end

    net.cs_vec = cs_vec
    net.cs_flag_vec = cs_vec2cs_flag_vec(net, cs_vec)
    net.cs_dict = Dict([cs.idx for cs in cs_vec] .=> cs_vec)
    # n_node = nv(net)
    #net.distance_data = spzeros(n_node, n_node)
    #net.grade_data = spzeros(n_node, n_node)
    add_junction!(net; min_nei=4, cs_net_flag=true)
    return net
end


"""
Get charging stations and convert them into nodes in a network.

region: the region of the charging station, could be us or eu (europe)
"""
function get_charge_station(net::Network, scenario::String; truck_flag::Bool = false, fuel::Bool = false, carbon_dataset = CambiumDataset(), continent::Symbol = :us)
    carbon_dict = get_carbon_dict(carbon_dataset)

    latlons = get_charge_station(continent; truck_flag = truck_flag, fuel=fuel)
    node_idx_vec = closest_node(net, latlons; max_com_flag = true)
    unique!(node_idx_vec)
    cs_vec::Vector{ChargeStation} = Vector{ChargeStation}()
    # @show length(node_idx_vec)

    region_box_dict = get_region_box_dict(carbon_dataset)
    for i in 1:length(node_idx_vec)
        idx = node_idx_vec[i]
        node = getnode(net, idx)
        (;lat, lon) = node
        # @debug "process cs $i/$(length(cs_vec)), ($lat, $lon)"
        reg = latlon2carbonregion(lat, lon, continent, carbon_dataset, region_box_dict)    
        if isempty(reg)
            continue
        end
        cs = ChargeStation(;
            lat=node.lat, lon=node.lon, idx=idx, region=reg, 
            avg_price = get_price_stat(reg, scenario, carbon_dict, mean),
            min_price = get_price_stat(reg, scenario, carbon_dict, minimum),
            max_price = get_price_stat(reg, scenario, carbon_dict, maximum),
            idx_r = idx,
        )
        push!(cs_vec, cs)
    end
    return cs_vec
end

"""
change the carbondataset of a cs_vec
"""
function change_carbon_dataset!(cs_vec::Vector, carbon_dataset::AbstractCarbonDataset, scenario, continent::Symbol = :us)
    carbon_dict = get_carbon_dict(carbon_dataset)
    region_box_dict = get_region_box_dict(carbon_dataset)
    for cs in cs_vec
        (;lat, lon) = cs
        reg = latlon2carbonregion(lat, lon, continent, carbon_dataset, region_box_dict)    
        cs.region = reg
        cs.avg_price = get_price_stat(reg, scenario, carbon_dict, mean)
        cs.min_price = get_price_stat(reg, scenario, carbon_dict, minimum)
        cs.max_price = get_price_stat(reg, scenario, carbon_dict, maximum)
    end
    return cs_vec
end

"""
get the locations of the charging station
"""
function get_charge_station(continent::Symbol; truck_flag::Bool = false, fuel::Bool = false)
    if fuel
        osmpath = joinpath(g_net_data_path, String(continent), "$continent-fuel.osm")
    else
        osmpath = joinpath(g_net_data_path, String(continent), "$continent-charging_station.osm")
    end
	if !isfile(osmpath)
		error("File $osmpath does not exist.")
		return  nothing
	end
    latlons = Vector{LatLon}()
    function istruck(xml_node)
        if truck_flag
            tag_dict = get_tag_dict(xml_node)
            return haskey(tag_dict, "truck") && (tag_dict["truck"] == "yes")
        else
            return true
        end
    end
    open(EzXML.StreamReader, osmpath) do reader
		for typ in reader
			if reader.name == "node" && typ == EzXML.READER_ELEMENT
				xml_node = expandtree(reader)
				## todo: elevation data missed in osm data
                if istruck(xml_node)
                    lat = Base.parse(Float64, xml_node["lat"])
                    lon = Base.parse(Float64, xml_node["lon"])
                    latlon = LatLon(lat, lon)
                    push!(latlons, latlon)
                end
			end
		end
	end # open the reader file

    return latlons
end

"""
Get average price for each region.
"""
function get_price_stat(region::AbstractString, scenario::String, carbon_dict, func::Function = mean, )
    # @show region scenario
    ta = carbon_dict[scenario][region]
    ta = TS.to( ta, DateTime(2025, 1, 1))
    pv = values(ta)
    idx = findall(x->(x!=0 && !isnan(x)), pv)
    return func(pv[idx])
end


"""
sparsify the charging stations such that they are not too close
node_vec: we also want the cs is far away from the nodes in node_vec
"""
function sparsify_cs_vec(net::Network, cs_vec::Vector{ChargeStation}, min_dis::Float64, node_vec::Vector{Int} = Int[])
    new_cs_vec::Vector{ChargeStation} = []

    for ics in 2:length(cs_vec)
        cs_to_add = cs_vec[ics] 
        
        u = cs_to_add.idx
        ## Do not an endpoint
        if is_endpoint(net, u)
            continue  
        end
        node_vec1 = vcat(node_vec, [cs.idx for cs in new_cs_vec])
        dis_vec = (distance3d(net, u, v) for v in node_vec1)
        if minimum(dis_vec; init=Inf) > min_dis
            # @debug "" u minimum(dis_vec)/1e3 min_dis/1e3
            push!(new_cs_vec, cs_to_add)
        end
    end

    return new_cs_vec 
end