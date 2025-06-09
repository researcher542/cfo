"""
It seems that 
"""

import E2PilotCore as ep
import E2PilotCFO as cfo
cfo.init()
using PyPlot
e2plt = Base.get_extension(cfo, :PlotExt)
e2plt.set_pyplot_style_no_latex()

if !isdefined(Main, :prob)
    prob = cfo.get_test_prob("eu")
end

(;min_lat, max_lat, min_lon, max_lon) = ep.get_eu_osm_box()

function plot_poly(ax, poly::ep.Polygon)
    coor = poly.coordinates[1]
    lon, lat = coor[1]
    if lon < min_lon || lon > max_lon || lat < min_lat || lat > max_lat
        return
    end
    poly_py = matplotlib.patches.Polygon(coor, closed=true, fill=false, edgecolor="black", linewidth=0.5)
    ax.add_patch(poly_py)
end

function plot_poly(ax, poly::ep.MultiPolygon)
    for p in poly.coordinates
        plot_poly(ax, p)
    end 
end

function plot_continent(ax, sym)
    region_data = ep.g_region_dict[sym]

    for feature in region_data.features
        plot_poly(ax, feature.geometry) 
    end
end


net = prob.net

fig, ax = plt.subplots()

function plot_od_data(ax)
    od_vec = cfo.read_od_data_vec(prob)
    x_vec = Float64[]
    y_vec = Float64[]
    for (i, od) in enumerate(od_vec)
        if !(i <= 200 && i > 100) 
            continue
        end
        nd1 = ep.getnode(net, od.src)
        nd2 = ep.getnode(net, od.des)
        dis = ep.distance2d(net, od.src, od.des) / ep.g_meter_per_mile
        push!(x_vec, nd1.lon)
        push!(y_vec, nd1.lat)
        # push!(x_vec, nd2.lon)
        # push!(y_vec, nd2.lat)
        ax.scatter(nd1.lon, nd1.lat, s=0.5, c="red", alpha=1.0)
        ax.scatter(nd2.lon, nd2.lat, s=0.5, c="blue", alpha=1.0)
        # if !( i <= 300 && i > 200) 
        
        println("i=$i Distance: ", dis)
    end

    # ax.scatter(x_vec, y_vec, s=0.5, c="red", alpha=1.0)
    @show length(x_vec) length(unique(x_vec))
    # global x_vec = x_vec
    
end

function plot_cs_vec(ax)
    cs_vec = prob.net.cs_vec
    x_vec = [cs.lon for cs in cs_vec]
    y_vec = [cs.lat for cs in cs_vec]
    ax.scatter(x_vec, y_vec, s=0.5, c="blue", alpha=1.0)
end



fig.set_size_inches(4, 5)
fig.tight_layout()
plot_continent(ax, :eu)
plot_od_data(ax)
# plot_cs_vec(ax)
e2plt.e2savefig(fig, "od_data")