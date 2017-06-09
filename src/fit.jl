module Optimize

include("orbits.jl")
include("observables.jl")

using Orbits: OrbitLibrary, BinnedLibrary
using Observables: KinematicData, binedges

function objective(lib::OrbitLibrary, data::KinematicData)
    binlib = BinnedLibrary(lib, data.edges)
    f = weights -> begin
        
    end
    return f
end

function fit(lib::OrbitLibrary, data::KinematicData)
    f = objective(lib, data)
    # optimization package tbd..., Ipopt or NLopt?
end

end
