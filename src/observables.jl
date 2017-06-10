module Observables

import Base: length, size, show, writedlm
using Formatting: format

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
end

function KinematicData(radii::Array{Float64, 1},
                       density::Array{Float64, 1},
                       unc_density::Array{Float64, 1},
                       dispersion::Array{Float64, 1},
                       unc_dispersion::Array{Float64, 1})
    nbins = length(radii)
    # get bin edges
    @assert nbins > 1
    edges = binedges(radii)
    # calculate mass in each bin
    dr = map(bin -> bin[2] - bin[1], edges)
    area = 2pi * radii .* dr
    mass = density .* area
    norm = sum(mass)
    mass = mass / norm
    unc_mass = unc_density .* area ./ norm
    KinematicData(radii, density, unc_density, dispersion, unc_dispersion, mass,
                  unc_mass, edges)
end

function KinematicData(filename::String)
    A = readdlm(filename)
    KinematicData(A[2], A[4], A[5], A[6], A[7])
end

show(io::IO, data::KinematicData) = println(format("KinematicData: {:d} bins", length(data)))
length(data::KinematicData) = length(data.radii)
size(data::KinematicData) = size(data.radii)
writedlm(filename::String, data::KinematicData) = begin
    rlow = [data.edges[i][1] for i in length(data)]
    rhigh = [data.edges[i][2] for i in length(data)]
    writedlm(filename, hcat(rlow, radii, rhigh, density, unc_density,
                            dispersion, unc_dispersion, mass, unc_mass))
end

function binedges(radii::Array{Float64, 1})
    nbins = size(radii, 1)
    edges = Array{Tuple{Float64, Float64}}(nbins)
    # first bin
    rhigh = (radii[1] + radii[2]) / 2    
    edges[1] = (max(0, 2 * radii[1] - rhigh), rhigh)
    # last bin
    rlow = (radii[end] + radii[end - 1]) / 2
    edges[end] = (rlow, 2 * radii[end] - rlow)
    for i in 2:(nbins - 1)
        edges[i] = ((radii[i - 1] + radii[i]) / 2, (radii[i] + radii[i + 1]) / 2)
    end
    return edges
end



end
