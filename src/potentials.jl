module Potentials

export PowerLawPotential, mass, density, potential

include("constants.jl")
using .Constants
c = Constants

abstract Potential
abstract SphericalPotential <: Potential

function density(phi::SphericalPotential, r::Real)
    error("Must use subtype of SphericalPotential with defined density!")
end

function mass(phi::SphericalPotential, r::Real)
    error("Numeric integration of density profile not implemented yet.")
end

function potential(phi::SphericalPotential, r::Real)
    M = mass(phi, r)
    return -c.G * M ./ r
end

immutable PowerLawPotential <: SphericalPotential
    #=
    r0 : scale radius
    rho0 : scale density
    alpha : power law index of density
    =#
    r0::Real
    rho0::Real
    alpha::Real
end

function mass(phi::PowerLawPotential, r::Real)
    r0 = phi.r0
    rho0 = phi.rho0
    alpha = phi.alpha
    4pi * rho0 * r0 ^ 3 / (3 + alpha) * (r / r0) .^ (3 + alpha)
end

function density(phi::PowerLawPotential, r::Real)
    r0 = phi.r0
    rho0 = phi.rho0
    alpha = phi.alpha
    rho0 * (r / r0) .^ alpha
end

end
