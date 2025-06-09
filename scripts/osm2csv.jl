"""
Convert the OSM data to CSV format.

node.csv has columns: name, lat, lon, ele, id

edge.csv has columns: src, des, name, speed_pos, speed_neg (in mph)
"""

import E2PilotCore as ep
import E2PilotCFO as cfo
import E2PilotServer as es
using CSV, HTTP, DataFrames, Geodesy, Statistics, Graphs, PyPlot
e2plt = Base.get_extension(cfo, :PlotExt)
e2plt.set_pyplot_style_no_latex()

function plotnet_bnd(net)
    fig, ax = plt.subplots()
    min_lat = Inf
    max_lat = -Inf
    min_lon = Inf
    max_lon = -Inf
    
    x_vec = Float64[]
    y_vec = Float64[]
    for nd in net.nodesdata
        lat = nd.lat
        lon = nd.lon
        if lat < max_lat && lat > min_lat && lon < max_lon && lon > min_lon
            ## if it is inside the bounding box
            # continue
        end
        push!(x_vec, lon)
        push!(y_vec, lat)
        
        min_lat = min(min_lat, lat)
        max_lat = max(max_lat, lat)
        min_lon = min(min_lon, lon)
        max_lon = max(max_lon, lon)
    end
    @show length(x_vec) length(y_vec)
    ax.scatter(x_vec, y_vec, s=1.0, c="red", alpha=0.5)
     
    fig.set_size_inches(12, 8)
    fig.tight_layout()
    e2plt.e2savefig(fig, "net_bnd.$(net.region)")
    close(fig)
end




function isvalid_ele(ele)
    if isnan(ele) || ele == 0.0 || abs(ele) > 10_000.0 || ele == -1.0
        return false
    end 
    return true
end

function request_elevation_df!(df)
    nrow_per_batch = 4000
    # ll = LatLon(48.8566, 2.3522) # Paris
    start_idx = findfirst(row -> !isvalid_ele(row.ele), eachrow(df))
    if isnothing(start_idx)
        @info "All nodes have been processed, no need to request elevation data."
        return
    end
    nrow = size(df, 1)
    end_idx = min(start_idx + nrow_per_batch, nrow)
    t_begin = time()
    t_last = time()
    while end_idx <= nrow + 1
        perc = round(start_idx / nrow * 100, digits=2)
        idx_vec = filter(idx -> !isvalid_ele(df.ele[idx]), start_idx:end_idx)
        @debug "processing $(start_idx)-th to $(end_idx)-th node, epsipled time $(time() - t_begin), progress: $(perc)%"
        ll_vec = [
            LatLon(df.lat[i], df.lon[i]) for i in idx_vec
        ]
        @debug "len(idx_vec) = $(length(idx_vec))"
        if !isempty(idx_vec)
            @time ele_vec = es.request_elevation(ll_vec, es.ElevationJl(); dataset="eudem")
            for i in 1:length(idx_vec)
                irow = idx_vec[i]
                ele = ele_vec[i]
                if isvalid_ele(ele)
                    df.ele[irow] = ele
                    # @debug "node $(irow) has valid elevation data $(ele)"
                else
                    df.ele[irow] = ele
                    @warn "node $(irow) has invalid elevation data $(ele)" # maxlog=10
                end
            end
        end

       
        # for i in start_idx:end_idx
        #     df.ele[i] = ele_vec[i - start_idx + 1]
        # end 
        if end_idx == nrow
            break
        end

        start_idx = end_idx + 1
        end_idx = min(start_idx + nrow_per_batch, nrow)
        t_eps = time() - t_last

        if t_eps > 60.0
            @info "Elapsed time: $(time() - t_begin) seconds, saving the temporary data to $(csv_data_path)"
            @time CSV.write(csv_data_path, df)
            GC.gc(true)
            t_last = time()
        end
        
        # break
    end

    @time CSV.write(csv_data_path, df)
    
end

ENV["JULIA_DEBUG"] = "Main,E2PilotCore"
# prob = cfo.get_test_prob("TXOH"); net = prob.net
region = ["oldmap", "eu", "eu.merged.v1", "eu.merged.v2"][4]
if !isdefined(Main, :net)
    # global net = ep.readmapdata(region; osm_flag=true);
    @time net = ep.readmapdata(region);
end
# ep.net2csv(net) ###############
# ep.net2nodecsv(net)
# ep.net2edgecsv(net)

### add the elevation data to the node data
csv_data_dir = joinpath(ep.k_data_path, "map", region)
csv_data_path = joinpath(csv_data_dir, "node.csv")
csv_edge_data_path = joinpath(csv_data_dir, "edge.csv")
if !isdefined(Main, :node_df)
    # @time node_df = DataFrame(CSV.File(csv_data_path))
    # @time edge_df = DataFrame(CSV.File(csv_edge_data_path))
end

function printnetstat(net)
    waydata_vec = collect(values(net.waydata))
    @info "Network statistics for region $(net.region)"
    for sym in [:distance, :grade]
        data_vec = [getproperty(way, sym) for way in waydata_vec]
        println("$(sym)" * "*"^80)
        @show sym
        @show mean(data_vec) minimum(data_vec) maximum(data_vec)
    end

    println("misc stat" * "*"^80)
    short_edge_vec = findall(w -> w.distance < 200.0, waydata_vec)
    @show length(short_edge_vec)
    onedegree_node_vec = findall( u -> (indegree(net, u) == 1 || outdegree(net, u) == 1), 1:nv(net))
    @show length(onedegree_node_vec)
    @show length(net.max_com)
    

    println("latlon stat" * "*"^80)
    lat_vec = (nd.lat for nd in net.nodesdata)
    lon_vec = (nd.lon for nd in net.nodesdata)
    @show minimum(lat_vec) maximum(lat_vec)
    @show minimum(lon_vec) maximum(lon_vec)
    @show ne(net) nv(net)
end


kwargs = Dict(
    :min_dis => 200.0,
    :grade_tol => 2.5e-3,
    :region=> "eu.merged.v2"
)
if !isdefined(Main, :net1) || true
end

@time net1 = ep.merge_one_degree_node(net; kwargs...)

# @time net2 = ep.induced_network_by_edge(net1)
# @time net3 = ep.induced_network_by_node(net1, 100)
printnetstat(net)
printnetstat(net1)

plotnet_bnd(net1)

# kwargs[:region] = "eu.merged.v3"
# @time net2 = ep.merge_one_degree_node(net1; kwargs...); printnetstat(net2)

# printnetstat(net2)

# net2 = ep.net_merge_nodes(net1, "eu.merged.v2", 1_000.0); printnetstat(net2)

# ep.net2csv(net3)

@info "done."
# net1 =ep.net_merge_nodes(net, region, 20_000.0)
# request_elevation_df!(node_df)

