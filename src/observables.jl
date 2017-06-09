module Observables

abstract Observable

type KinematicData <: Observable
    radii::Array{Float64, 1}
    density::Array{Float64, 1}
    unc_density::Array{Float64, 1}
    dispersion::Array{Float64, 1}
    unc_dispersion::Array{Float64, 1}
    mass::Array{Float64, 1}
    unc_mass::Array{Float64, 1}
    edges::Array{Tuple{Float64, Float64}}
    function KinematicData(radii::Array{Float64, 1},
                           density::Array{Float64, 1},
                           unc_density::Array{Float64, 1},
                           dispersion::Array{Float64, 1},
                           unc_dispersion::Array{Float64, 1})
        # get bin edges
        @assert nbins > 1
        edges = binedges(radii)
        # calculate mass in each bin
        dr = map(bin -> edges[2] - edges[1], edges)
        area = 2pi * radii .* dr
        mass = density .* area
        unc_mass = unc_density .* area
        new(radii, density, unc_density, dispersion, unc_dispersion, mass, unc_mass, edges)
    end
end

function KinematicData(file::String)
    KinematicData(readdlm(file)...)
end

function binedges(radii::Array{Float64, 1})
    nbins = size(radii, 1)
    edges = Array{Tuple{Float63, Float64}}(nbins)
    # first bin
    rhigh = (radii[1] + radii[2]) / 2
    edges[1] = (2 * radii[1] - rhigh, rhigh)
    # last bin
    rlow = (radii[end] + radii[end - 1]) / 2
    edges[end] = (rlow, 2 * radii[end] - rlow)
    for i in 2:(nbins - 1)
        edges[i] = ((radii[i - 1] + radii[i]) / 2, (radii[i] + radii[i + 1]) / 2)
    end
    return edges
end

import Base.size
size(data::KinematicData) = size(data.radii)

end
