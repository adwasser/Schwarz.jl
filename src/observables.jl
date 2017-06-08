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
    bins::Array{Tuple{Float64, Float64}}
    function KinematicData(radii::Array{Float64, 1},
                           density::Array{Float64, 1},
                           unc_density::Array{Float64, 1},
                           dispersion::Array{Float64, 1},
                           unc_dispersion::Array{Float64, 1})
        # get bin edges
        nbins = size(radii)[1]
        @assert nbins > 1
        bins = Array{Tuple{Float64, Float64}}(nbins)
        # first bin
        rhigh = (radii[1] + radii[2]) / 2
        bins[1] = (2 * radii[1] - rhigh, rhigh)
        # last bin
        rlow = (radii[end] + radii[end - 1]) / 2
        bins[end] = (rlow, 2 * radii[end] - rlow)
        for i in 2:(nbins - 1)
            bins[i] = ((radii[i - 1] + radii[i]) / 2, (radii[i] + radii[i + 1]) / 2)
        end
        
        # calculate mass in each bin
        dr = map(bin -> bin[2] - bin[1], bins)
        area = 2pi * radii .* dr
        mass = density .* area
        unc_mass = unc_density .* area
        new(radii, density, unc_density, dispersion, unc_dispersion, mass, unc_mass, bins)
    end
end

function KinematicData(file::String)
    KinematicData(readdlm(file)...)
end

import Base.size
size(data::KinematicData) = size(data.radii)

end
