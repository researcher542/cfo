"""
self-implemented GeoJSON
"""

@with_kw struct Feature{T,S}
    geometry::T
    properties::S
end

@with_kw struct FeatureCollection{T}
    features::Vector{T}
end

@with_kw struct Polygon{T}
    coordinates::Vector{T}
end

@with_kw struct MultiPolygon{T}
    coordinates::Vector{T}
end

"""
convert json object to a Polygon or MultiPolygon
"""
function json2geometry(obj::JSON3.Object)
    obj_type = obj[:type]
    if obj_type == "Feature"
        geo = json2geometry(obj[:geometry])
        return Feature(geo, obj[:properties])
    elseif obj_type == "Polygon"
        coor = Base.copy(obj[:coordinates])
        coor1 = coor
        # coor1 = [
        #     [SA[xy[1], xy[2]] for xy in cor]
        #     for cor in coor
        # ]
        return Polygon(coor1)
    elseif obj_type == "MultiPolygon"
        coor_vec = Base.copy(obj[:coordinates])
        poly_vec = [Polygon(co) for co in coor_vec]
        return MultiPolygon(poly_vec)
    else
        error("Unknown type $obj_type")
        return nothing
    end
end

function load_region(region::Symbol)
    data_path = joinpath(k_root_dir,"data","region","$region.geojson")
    f = read(data_path)
    json = JSON3.read(f)
    # json1 = Base.copy(json)
    @assert json.type == "FeatureCollection"
    feature_vec = [json2geometry(feature) for feature in json.features]
    return FeatureCollection(feature_vec)

    # for feature in json1[:features]
    #     feature[:geometry] = json2geometry(feature[:geometry])
    # end
    # object = Base.copy(js)
    # object_type = get(object, :type, nothing)
    # @assert object_type == "FeatureCollection"

    # data = GeoJSON.FeatureCollection(DictWapper(object))
    # copy json3 to ensure type stability.
    # data = GeoJSON.read(f)
    # return json1
end