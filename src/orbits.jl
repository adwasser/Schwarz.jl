module Orbits

include("potentials.jl")
import .Potentials: mass, density, potential, SphericalPotential

type Orbit
    #=
    e : specific energy
    j : specific angular momentum
    =#
    e::Real
    j::Real
end

type OrbitLibrary
    phi::SphericalPotential
end

end
