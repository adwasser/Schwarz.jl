module Orbits

using Roots.fzero
using DifferentialEquations: ODEProblem, solve, Vern6

using Schwarz.Constants
c = Constants
using Schwarz.Potentials: mass, density, potential, SphericalPotential

type Orbit
    #=
    Orbit in spherical potential, restricted to x-y plane

    e : specific energy
    j : specific angular momentum
    phi : SphericalPotential
    =#
    e::Real
    j::Real
    phi::SphericalPotential
    x::Array
    y::Array
    vx::Array
    vy::Array
    t::Array

    Orbit(w0, phi) = begin
        pos = w0[1:2]
        vel = w0[3:4]
        r = norm(pos)
        v = norm(vel)
        j = norm(cross(vcat(pos, [0]), vcat(vel, [0])))
        e = potential(phi, r) + v ^ 2 / 2        
        g(t, w) = begin
            n = Int(size(w)[1] / 2)
            pos = w[1:n]
            vel = w[n + 1:end]
            r = sqrt(sum(pos .^ 2))
            vcat(vel, (-c.G * mass(phi, r) / r ^ 3) .* pos)
        end
        # integrate for 10 Gyr
        t_int = 10 * c.Gyr
        dt = 1e-2 * c.Gyr
        prob = ODEProblem(g, w0, (0, t_int))
        sol = solve(prob, Vern6())
        t = collect(0:dt:t_int)
        w = sol(t)
        x = map(w -> w[1], w)
        y = map(w -> w[2], w)
        vx = map(w -> w[3], w)
        vy = map(w -> w[4], w)
        new(e, j, phi, x, y, vx, vy, t)        
    end
    
    Orbit(e, j, phi) = begin
        x0 = fzero(r -> e - potential(phi, r) - j ^ 2 / r ^ 2 / 2,
                   1e-3 * c.kpc, 1e3 * c.kpc)
        vy0 = sqrt(2 * (e - potential(phi, r0)))
        w0 = [x0, 0, 0, vy0]
        Orbit(w0, phi)
    end
end

type OrbitLibrary
    phi::SphericalPotential
end


end
