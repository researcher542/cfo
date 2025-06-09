
function readnodes_csv!(net::Network, region::String = "test")
	csvfile = region2nodecsv(region)	
	df = DataFrame(CSV.File(csvfile))
	idx = 0
	for row in eachrow(df)
		idx += 1
		id = row.id
		node = Node(row.lat, row.lon, row.elevation, row.id, idx)
		addnode!(net, node)
	end
	return net
end

"""
@Note: in osm data, a junction node is not a actual junction... Should look the trunk_link in the way representation.
"""
function readnodes_osm!(net::Network, region::String = "test")
	osmpath = region2osmfile(region)
	if !isfile(osmpath)
		@error "File $osmpath does not exist."
		return missing
	end
    @info "reading nodes from osm file $(osmpath)"
	open(EzXML.StreamReader, osmpath) do reader
		idx = 0
		for typ in reader
			#@show typ, reader.name, reader.content, reader.depth
			if reader.name == "node" && typ == EzXML.READER_ELEMENT
				#@show typ, reader.name, reader.content, reader.depth
				xml_node = expandtree(reader)
				idx += 1
				if idx % 100000 == 0
					@debug "processing $(idx)-th node"
				end
				ele = NaN
				## todo: elevation data missed in osm data
				id = str2int(xml_node["id"])
				tag_dict = get_tag_dict(xml_node)
				net_node = Node(xml_node["lat"], xml_node["lon"], ele, id, idx)
				addnode!(net, net_node)
				#add_vertex!(net.graph)	
				#push!(net.nodesdata, net_node)
				#net.nodesdata[idx] = net_node
				#@debug tag_dict
				#set_props!(net.graph, tag_dict)
				for (k,v) in tag_dict
					#set_prop!(net.graph, idx, Symbol(k), v)
					#if occursin("junction", v)
					#	push!(net.junctions, idx)
					#end
				end
				#net.id2idx[id] = idx
			end
		end
	end # open the reader file
	return net
end