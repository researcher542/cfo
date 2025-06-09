"""
Provide some geometry functions
geojson data of us is from https://eric.clst.org/tech/usgeojson/
geojson data of eu is from 
    - (used) https://ec.europa.eu/eurostat/web/gisco/geodata/administrative-units/countries CNTR_RG_60M_2024_4326
    - https://github.com/leakyMirror/map-of-europe/blob/master/GeoJSON/europe.geojson

"""

# import GeoJSON
using PolygonOps
using Geodesy
using GeoInterface 

function get_eu_osm_box()
    g_eu_osm_box = (min_lat = 36.06, max_lat = 61.88, min_lon = -9.4, max_lon = 23.09)
    return g_eu_osm_box
end


function inpolygon(lat::Real, lon::Real, multi_poly::MultiPolygon)
    for poly in multi_poly.coordinates
        flag = inpolygon(lat, lon, poly)
        if flag != 0
            return flag
        end
    end
    return 0
end

function inpolygon(lat::Real, lon::Real, poly::Polygon)
    # return inpolygon_fast(lat, lon, Base.copy(poly.coordinates[1]))
    # coordinates = Base.copy(poly.coordinates[1])
    return inpolygon_fast(lat, lon, poly.coordinates[1])
end

"""
self implementation of the inpolygon, with fast check
"""
function inpolygon_fast(lat::Real, lon::Real, poly)
    x = lon
    y = lat
    # @show typeof(poly)
    xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
    @inbounds for ixy in 1:length(poly)
        # xy::JSON3.Array{Float64, Vector{UInt8}, SubArray{UInt64, 1, Vector{UInt64}, Tuple{UnitRange{Int64}}, true}} = poly[ixy]
        xy = poly[ixy]
        # xy = poly[ixy]
        @inbounds xx::Float64 = xy[1]
        @inbounds yy::Float64 = xy[2]
        xmin::Float64 = min(xmin, xx)
        xmax::Float64 = max(xmax, xx)
        ymin::Float64 = min(ymin, yy)
        ymax::Float64 = max(ymax, yy)
    end
    # xmin,xmax = minimum(x_vec), maximum(x_vec)
    # ymin,ymax = minimum(y_vec), maximum(y_vec)
    if (x < xmin) || (x > xmax) || (y < ymin) || (y > ymax)
        return 0
    end
    # poly1 = [poly[1:2:end]..., poly[end]]
    flag = PolygonOps.inpolygon((lon, lat), poly)
    return flag
end

global g_region_dict = nothing
function init_geo()
    global g_regions = [:china, :us, :eu]
    if isnothing(g_region_dict)
        global g_region_dict = Dict(g_regions .=> load_region.(g_regions))
    end
    global g_hk_poly = subregion_poly(:china, "香港特别行政区")
end

function getname(region, poly)
    if region == :china
        return poly.properties["name"]
    elseif region == :us
        return poly.properties["name"]
    elseif region == :eu
        return poly.properties["NAME_ENGL"]
    else
        @warn "region $region not supported"
        # @assert false
    end
    
end

"""
Get a polygon of a subregion
"""
function subregion_poly(region::Symbol, subregion::String)
    if !(region in keys(g_region_dict))
        @warn "$region not in the region_dict"
        return nothing
    end
    region_data = g_region_dict[region]
    idx_vec = findall( 
        x-> (
            getname(region, x) == subregion
            ), 
        region_data.features
        )
    if isempty(idx_vec)
        @warn "subregion $subregion not found in $region"
        @assert false
        return nothing
    end
    idx = idx_vec[1]
    # return region_data.geometry[idx]
    return region_data.features[idx].geometry
    # return region_data[:features][idx][:geometry]
end


function isinregion(lat::Real, lon::Real, region::Symbol)
    return !(isnothing(latlon2subregion(lat, lon, region)))
end

function latlon2subregion(lat::Real, lon::Real, region::Symbol, verbose::Int = 0)
    if !(region in keys(g_region_dict))
        @warn "$region not in the region_dicts"
        return nothing
    end
    region_data = g_region_dict[region]
    for (ipoly, feature) in enumerate(region_data.features)
        poly = feature.geometry
        # name = poly.properties["name"]
        # geometry = poly.geometry
        flag = inpolygon(lat, lon, poly)
        if flag != 0
            # p = feature.properties
            name = getname(region, feature)
            return name
            # return feature.properties.name
            # return region_data.name[ipoly]
            # return name
        end
    end
    if verbose > 0
        @warn "($lat, $lon) not in the region $region"
    end
    return nothing
end

isinregion(p::LatLon, region::Symbol) = isinregion(p.lat, p.lon, region)

function outofchina_fast(lat::Real, lon::Real)
    if (lon < 72.004) | (lon > 137.8347) 
        return true
    end
    if (lat < 0.8293) | (lat > 55.8271) 
        return true
    end
    return false
end

function isinchina(lat::Real, lon::Real)  
    if outofchina_fast(lat, lon)
        return false
    end
    return isinregion(lat, lon, :china)
end
isinchina(p::LatLon) = isinchina(p.lat, p.lon)

function isinHK(lat::Real, lon::Real)
    flag = inpolygon(lat, lon, g_hk_poly)
    return flag != 0
end
isinHK(latlon::LatLon) = isinHK(latlon.lat, latlon.lon)