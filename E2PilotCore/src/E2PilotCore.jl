module E2PilotCore


# requirement by the network
using Parameters, JSON3, TimerOutputs, Graphs, DataStructures, Infiltrator, CSV, DataFrames, Roots, FastPriorityQueues, ProgressMeter, StaticArrays, Dates, Dictionaries, Polyester, ZipArchives
import Geodesy: LatLon, LLA, euclidean_distance
import Distances: Haversine
import Distances

import NearestNeighbors as NN

# requirement by the planning
using EzXML, Interpolations


include("const.jl")
include("utils/E2Utils.jl")
include("network/E2Network.jl")
include("route/route.jl")
include("plan/E2Plan.jl")
include("oddata/oddata.jl")

# greet() = print("Hello World!")
function init()
    @info "Initilizing E2PilotCore"
    init_geo()
end

# include("precompile.jl")

end # module E2PilotCore
