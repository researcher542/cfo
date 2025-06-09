"""
Convert the net to CSV format.

node.csv has columns: name, lat, lon, ele, id

edge.csv has columns: src, des, name, speed_pos, speed_neg (in mph)
"""
function net2csv(net::AbsNet)
    net2nodecsv(net)
    net2edgecsv(net)
end

function net2edgecsv(net::AbsNet)
    df_edge_dict = OrderedDict()
    waydata_vec = collect(values(net.waydata))
    for sym0 in [:src, :des]
        sym = sym0
        df_edge_dict[sym0] = [getproperty(w, sym) for w in waydata_vec]
    end

    default_spd_mph = 60.0
    df_edge_dict[:name] = ["" for _ in 1:length(df_edge_dict[:src])]
    df_edge_dict[:speed_pos] = [default_spd_mph for _ in 1:length(df_edge_dict[:src])]
    df_edge_dict[:speed_neg] = [default_spd_mph for _ in 1:length(df_edge_dict[:src])]

    # filepath = region2edgecsv(net.region)
    @info "Writing edge data to CSV file..."
	datadir = joinpath(k_root_dir, "data", "map", net.region)
    filepath = joinpath(datadir, "edge.csv")
    df = DataFrame(df_edge_dict)
    CSV.write(filepath, df)
end

function net2nodecsv(net::AbsNet)
    df_node_dict = OrderedDict()
    
    for sym in [:name, :lat, :lon, :ele, :id]
        df_node_dict[sym] = [getproperty(nd, sym) for nd in net.nodesdata]
    end

    @info "Writing node data to CSV file..."
	datadir = joinpath(k_root_dir, "data", "map", net.region)
    if isdir(datadir) == false
        mkpath(datadir)
    end


    filepath = joinpath(datadir, "node.csv")
    df = DataFrame(df_node_dict)
    CSV.write(filepath, df)
    
end