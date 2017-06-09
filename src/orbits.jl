module Orbits

using Roots.fzero

# using DifferentialEquations: ODEProblem, solve, AbstractODEAlgorithm
using DifferentialEquations
default_solver = DifferentialEquations.Vern9()

using Schwarz.Constants
c = Constants
using Schwarz.Potentials: mass, density, potential, SphericalPotential
using Schwarz.Rotations: Rx, Ry, Rz

type Orbit
    #=
    Orbit in spherical potential, restricted to x-y plane

    e : specific energy
    j : specific angular momentum
    i : inclination (rad)
    omega : longitude of ascending node (rad)
    phi : SphericalPotential
    =#
    e::Real
    j::Real
    i::Real
    omega::Real
    phi::SphericalPotential
    x::Array
    y::Array
    z::Array
    vx::Array
    vy::Array
    vz::Array
    t::Array
end

function Orbit(w0::Array, i::Real, omega::Real,
                   phi::SphericalPotential, alg::AbstractODEAlgorithm)
    #=
    w0 : intial condition, [x, y, vx, vy]
    i : inclination (rad)
    omega : longitude of ascending node (rad)
    phi : SphericalPotential
    =#
    pos = w0[1:2]
    vel = w0[3:4]
    r = norm(pos)
    v = norm(vel)
    j = norm(cross(vcat(pos, [0]), vcat(vel, [0])))
    e = potential(phi, r) + v ^ 2 / 2        
    g(t, w, dw) = begin
        n = 2
        pos = w[1:n]
        vel = w[(n + 1):end]
        r = sqrt(sum(pos .^ 2))
        dw .= vcat(vel, (-c.G * mass(phi, r) / r ^ 3) .* pos)
    end
    # integrate for 10 Gyr
    t_int = 10 * c.Gyr
    dt = 1e-3 * c.Gyr
    prob = ODEProblem(g, w0, (0, t_int))
    sol = solve(prob, alg)
    t = collect(0:dt:t_int)
    w = sol(t)
    x = map(w -> w[1], w)
    y = map(w -> w[2], w)
    z = map(w -> 0.0, w)
    vx = map(w -> w[3], w)
    vy = map(w -> w[4], w)
    vz = map(w -> 0.0, w)
    # rotate with i and omega
    pos = Rz(omega) * Rx(i) * transpose([x y z])
    x = pos[1, :]
    y = pos[2, :]
    z = pos[3, :]
    vel = Rz(omega) * Rx(i) * transpose([vx vy vz])
    vx = vel[1, :]
    vy = vel[2, :]
    vz = vel[3, :]
    Orbit(e, j, i, omega, phi, x, y, z, vx, vy, vz, t)        
end

Orbit(w0::Array, i::Real, omega::Real, phi::SphericalPotential) = Orbit(w0, i, omega, phi, default_solver)

Orbit(w0::Array, phi::SphericalPotential, alg::AbstractODEAlgorithm) = Orbit(w0, 0, 0, phi, alg)

Orbit(w0::Array, phi::SphericalPotential) = Orbit(w0, phi, default_solver)
    
function rotate(orbit::Orbit, i::Real, omega::Real)
    pos = transpose([orbit.x orbit.y orbit.z])
    vel = transpose([orbit.vx orbit.vy orbit.vz])
    pos = Rz(omega) * Rx(i) * pos
    vel = Rz(omega) * Rx(i) * vel
    x, y, z = pos[1, :], pos[2, :], pos[3, :]
    vx, vy, vz = vel[1, :], vel[2, :], vel[3, :]
    
    return Orbit(orbit.e, orbit.j, orbit.i + i, orbit.omega + omega, orbit.phi,
                 x, y, z, vx, vy, vz, orbit.t)
end

function rotate!(orbit::Orbit, i::Real, omega::Real)
    pos = transpose([orbit.x orbit.y orbit.z])
    vel = transpose([orbit.vx orbit.vy orbit.vz])
    pos = Rz(omega) * Rx(i) * pos
    vel = Rz(omega) * Rx(i) * vel
    x, y, z = pos[1, :], pos[2, :], pos[3, :]
    vx, vy, vz = vel[1, :], vel[2, :], vel[3, :]
    orbit.x = x
    orbit.y = y
    orbit.z = z
    orbit.i += i
    orbit.omega += omega
    return nothing
end

import Base: size, writedlm
size(orbit::Orbit) = size(orbit.t)
writedlm(filename::String, o::Orbit) = begin
    writedlm(filename, hcat(o.t / c.Gyr, o.x / c.kpc, o.y / c.kpc, o.z / c.kpc,
                            o.vx / c.km_per_s, o.vy / c.km_per_s,
                            o.vz / c.km_per_s))
end

function integrals(orbit::Orbit)
    # Calculate specific energy and specific angular momentum
    pos = [orbit.x orbit.y orbit.z]
    vel = [orbit.vx orbit.vy orbit.vz]
    r = map(i -> norm(pos[i, :]), 1:size(pos, 1))
    v = map(i -> norm(vel[i, :]), 1:size(vel, 1))
    e = map(i -> potential(orbit.phi, r[i]) + v[i]^2 / 2, 1:size(pos, 1))
    j = map(i -> norm(cross(pos[i, :], vel[i, :])), 1:size(pos, 1))
    return e, j
end

type OrbitLibrary
    phi::SphericalPotential
    orbits::Array{Orbit}
    function OrbitLibrary(Nr::Int, Nv::Int, Ni::Int, Nomega::Int,
                          rmax::Real, phi::SphericalPotential)
        dr = rmax / Nr
        radii = linspace(dr, rmax, Nr)
        orbits = Array{Orbit}(Nr * Nv * Ni * Nomega)
        count = 1
        for r in radii
            vmax = 0.9 * sqrt(-2 * potential(phi, r))
            dv = vmax / Nv
            velocities = linspace(dv, vmax, Nv)
            for v in velocities
                w0 = [r, 0, 0, v]
                orbit = Orbit(w0, phi)
                di = pi / 2.0 / Ni
                inclinations = linspace(di, pi / 2.0, Ni)
                for i in inclinations
                    domega = 2pi / Nomega
                    longitudes = linspace(domega, 2pi, Nomega)
                    for omega in longitudes
                        orbits[count] = rotate(orbit, i, omega)
                        count += 1
                    end
                end
            end
        end
        new(phi, orbits)
    end
end


end
