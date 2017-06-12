module Orbits

import Base: size, writedlm, show, length
using Formatting.format
# using DifferentialEquations
# default_solver = DifferentialEquations.Vern9()
using ODE

using Schwarz.Constants
c = Constants
using Schwarz.Potentials: mass, density, potential, SphericalPotential
using Schwarz.Rotations: Rx, Ry, Rz

type Orbit
    #=
    Orbit in spherical potential, restricted to x-y plane
    
    x, y : positions
    vx, vy : velocities
    =#
    # e::Real
    # j::Real
    # i::Real
    # omega::Real
    # phi::SphericalPotential
    t::Array
    x::Array
    y::Array
    vx::Array
    vy::Array
end

function Orbit(w0::Array, phi::SphericalPotential;
               t_int = 10 * c.Gyr, n_t = 4e3)
    #=
    Orbit in spherical potential, stored 2d coords

    w0 : intial condition, [x, y, vx, vy]
    phi : SphericalPotential
    =#
    g(t, w) = begin
        pos = w[1:2]
        vel = w[3:4]
        r = sqrt(sum(pos .^ 2))
        vcat(vel, (-c.G * mass(phi, r) / r ^ 3) .* pos)
    end
    dt = t_int / n_t
    t = collect(0:dt:t_int)
    tout, wout = ode78(g, w0, t; points=:specified)
    x = map(w -> w[1], wout)
    y = map(w -> w[2], wout)
    vx = map(w -> w[3], wout)
    vy = map(w -> w[4], wout)
    Orbit(t, x, y, vx, vy)        
end

Orbit(filename::String) = begin
    array = readdlm(filename)
    Orbit(array[:, 1] * c.Gyr, array[:, 2] * c.kpc, array[:, 3] * c.kpc,
          array[:, 4] * c.km_per_s, array[:, 5] * c.km_per_s)
end
show(io::IO, orbit::Orbit) = begin
    print(format("Orbit: x0 = {:.2f} kpc, v0 = {:.2f} km/s, t_int = {:.2f} Gyr",
                 orbit.x[1] / c.kpc, orbit.vy[1] / c.km_per_s, orbit.t[end] / c.Gyr))
end
size(orbit::Orbit) = size(orbit.x)
length(orbit::Orbit) = length(orbit.x)
writedlm(filename::String, o::Orbit) = begin
    writedlm(filename, hcat(o.t / c.Gyr, o.x / c.kpc, o.y / c.kpc,
                            o.vx / c.km_per_s, o.vy / c.km_per_s))
end

function rotate(orbit::Orbit, i::Real, omega::Real)
    rot = Rz(omega) * Rx(i)
    pos = rot * transpose(hcat(orbit.x, orbit.y, zeros(length(orbit))))
    vel = rot * transpose(hcat(orbit.vx, orbit.vy, zeros(length(orbit))))
    x, y, z = pos[1, :], pos[2, :], pos[3, :]
    vx, vy, vz = vel[1, :], vel[2, :], vel[3 ,:]
    return x, y, z, vx, vy, vz
end

function integrals(orbit::Orbit)
    # Calculate specific energy and specific angular momentum
    pos = [orbit.x orbit.y]
    vel = [orbit.vx orbit.vy]
    r = map(i -> norm(pos[i, :]), 1:length(pos))
    v = map(i -> norm(vel[i, :]), 1:length(vel))
    e = map(i -> potential(orbit.phi, r[i]) + v[i]^2 / 2, 1:length(pos))
    j = map(i -> norm(cross(pos[i, :], vel[i, :])), 1:length(pos))
    return e, j
end

"""
BinnedOrbit
-----------
mass : fraction mass in bin
dispersion : variance of line-of-sight velocity
index : into the orbit library for the orbit used to generate this binned orbit
"""
type BinnedOrbit
    mass::Array{Float64, 1}
    dispersion::Array{Float64, 1}
    index::Int
end

function BinnedOrbit(orbit::Orbit, i::Real, omega::Real,
                     edges::Array{Tuple{Float64, Float64}};
                     projection = [0.0, 0.0, 1.0])
    #=
    orbit : Orbit instance
    i : inclination (rad)
    omega : longitude of ascending node (rad)
    edges : lower and upper radii edges of bins
    projection : [x, y, z] of line-of-sight vector
    =#
    @assert size(projection, 1) == 3
    projection = projection / norm(projection)
    radii = map(edges -> (edges[1] + edges[2]) / 2, edges)
    mass = Array{Float64, 1}(size(edges, 1))
    dispersion = Array{Float64, 1}(size(edges, 1))
    x, y, z, vx, vy, vz = rotate(orbit, i, omega)
    pos = hcat(x, y, z)
    vel = hcat(vx, vy, vz)
    r = map(k -> norm(pos[k, :] - dot(pos[k, :], projection) * projection), 1:length(orbit))
    vlos = map(k -> dot(vel[k, :], projection), 1:length(orbit))
    for (i, (rlow, rhigh)) in enumerate(edges)
        inbin = (r .< rhigh) & (r .>= rlow)
        dispersion[i] = sum(vlos[inbin] .^ 2)
        # fraction of time spent in bin
        mass[i] = count(j -> j, inbin) / length(r)
    end
    BinnedOrbit(mass, dispersion, -1)
end

show(io::IO, orbit::BinnedOrbit) = println(format("BinnedOrbit: {:d} bins", length(orbit)))
size(orbit::BinnedOrbit) = size(orbit.mass)
length(orbit::BinnedOrbit) = length(orbit.mass)

type OrbitLibrary
    phi::SphericalPotential
    orbits::Array{Orbit}
end

function OrbitLibrary(Nr::Int, Nv::Int, rmax::Real, phi::SphericalPotential;
                      multithreaded = false)
    dr = rmax / Nr
    radii = logspace(log10(dr), log10(rmax), Nr)
    orbits = Array{Orbit}(Nr * Nv)
    # bins = Array{BinnedOrbit}(Nr * Nv * Ni * Nomega)
    initial_points = []
    for r in radii
        vmax = 0.9 * sqrt(-2 * potential(phi, r))
        dv = vmax / Nv
        velocities = linspace(dv, vmax, Nv)
        for v in velocities
            push!(initial_points, [r, 0, 0, v])
        end
    end
    if multithreaded
        f = pmap
    else
        f = map
    end 
    orbits = f(w0 -> Orbit(w0, phi), initial_points)
    OrbitLibrary(phi, orbits)
end

size(lib::OrbitLibrary) = size(lib.orbits)
length(lib::OrbitLibrary) = length(lib.orbits)
show(io::IO, lib::OrbitLibrary) = println(format("OrbitLibrary: {:d} orbits", size(lib.orbits)[1]))

type BinnedLibrary
    orbits::Array{BinnedOrbit}
    # unbinnedLib::OrbitLibrary
end

function BinnedLibrary(lib::OrbitLibrary, edges::Array{Tuple{Float64, Float64}},
                       Ni::Int, Nomega::Int; multithreaded = false)
    di = pi / 2.0 / Ni
    inclinations = linspace(di, pi / 2.0, Ni)
    domega = 2pi / Nomega
    longitudes = linspace(domega, 2pi, Nomega)
    inputs = [(lib.orbits[i], inclinations[j], longitudes[k]) for
              i in 1:length(lib),
              j in 1:length(inclinations),
              k in 1:length(longitudes)]
    if multithreaded
        f = pmap
    else
        f = map
    end
    orbits = f(i -> BinnedOrbit(inputs[i]..., edges), 1:length(inputs))
    # binlib = BinnedLibrary(orbits, lib)
    binlib = BinnedLibrary(orbits)
    # normalize mass
    for o in binlib.orbits
        o.mass = o.mass / length(binlib)
    end
    return binlib
end

function BinnedLibrary(Nr::Int, Nv::Int, Ni::Int, Nomega::Int,
                       rmax::Real, phi::SphericalPotential,
                       edges::Array{Tuple{Float64, Float64}})
    lib = OrbitLibrary(Nr, Nv, rmax, phi)
    BinnedLibrary(lib, edges, Ni, Nomega)
end

size(lib::BinnedLibrary) = size(lib.orbits)
length(lib::BinnedLibrary) = length(lib.orbits)
show(io::IO, lib::BinnedLibrary) = println(format("BinnedLibrary: {:d} orbits", size(lib.orbits)[1]))

end
