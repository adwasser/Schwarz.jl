module Optimize

include("orbits.jl")
include("observables.jl")

using Schwarz.Orbits: OrbitLibrary, BinnedLibrary
using Schwarz.Observables: KinematicData, binedges

function objective(lib::BinnedLibrary, data::KinematicData)
    #=
    From an orbit library and kinematic data, return the objective function,
    its gradient, constraints, etc.
    =#
    # lib = BinnedLibrary(lib, data.edges)
    # m = fractional mass times sigma^2
    # m_lib is n_orbit rows by n_bins columns matrix
    m_lib = [lib.orbits[i].mass[j] * lib.orbits[i].dispersion[j]
             for i in 1:length(lib), j in 1:length(data)]
    mass_lib = [lib.orbits[i].mass[j]
                for i in 1:length(lib), j in 1:length(data)]
    m_data = [data.mass[j] * data.dispersion[j] for j in 1:length(data)]
    dm_data = [sqrt((data.unc_mass[j] * data.dispersion[j])^2
                    + (data.unc_dispersion[j] * data.mass[j])^2)
               for j in 1:length(data)]
    u_data = 1 ./ dm_data .^ 2
    Q = Array{Float64, 2}(length(lib), length(lib))
    for i in 1:length(lib), j in 1:i
        if i == j
            Q[i, j] = 2 * sum(u_data .* m_lib[i, :].^2)
        else
            x = 4 * sum(u_data .* m_lib[i, :] .* m_lib[j, :])
            Q[i, j] = x
            Q[j, i] = x
        end
    end
    c = -2 * [sum(u_data .* m_data .* m_lib[i, :]) for i in 1:length(lib)]
    
    # objective
    function eval_f(weights::Vector{Float64})
        0.5 * dot(weights, Q * weights) + dot(c, weights)
    end
    
    # gradient of objective
    function eval_grad_f(weights::Vector{Float64}, grad_f::Vector{Float64})
        grad_f[:] = (Q * weights) .+ c
    end
    
    # constraints
    function eval_g(weights::Vector{Float64}, g::Vector{Float64})
        weighted_mass_lib = [dot(weights, mass_lib[:, j]) for j in 1:length(data)]
        # upper bound of 1
        g[1:(end - 1)] = abs(weighted_mass_lib .- data.mass) ./ data.unc_mass
        # upper and lower bound of 1
        g[end] = sum(weights)
    end

    # jacobian of constraints
    function eval_jac_g(weights::Vector{Float64}, mode,
                        rows::Vector{Int32}, cols::Vector{Int32},
                        values::Vector{Float64})
        #=
        rows: n_bins + 1 constraints, cols: n_orbit variables
        =#
        weighted_mass_lib = [dot(weights, mass_lib[:, j]) for j in 1:length(data)]
        sgn = map(j -> sign(data.mass[j] - weighted_mass_lib[j]),
                  1:length(data))
        if mode == :Structure
            idx = 1
            for row in 1:(length(data) + 1)
                for col in 1:length(lib)
                    rows[idx] = row
                    cols[idx] = col
                    idx += 1
                end
            end
        else
            count = 1
            for i in 1:length(lib)
                values[count:(count - 1 + length(data))] =
                    [sgn[j] * mass_lib[i, j] / data.unc_mass[j]
                     for j in 1:length(data)]
                values[count + length(data)] = 1.0
                count += length(data) + 1
            end
        end
    end

    function eval_h(weights::Vector{Float64}, mode,
                    rows::Vector{Int32}, cols::Vector{Int32},
                    obj_factor::Float64, lambda::Vector{Float64},
                    values::Vector{Float64})
        if mode == :Structure
            idx = 1
            for row in 1:(length(lib))
                for col in 1:row
                    rows[idx] = row
                    cols[idx] += 1
                end
            end
        else
            idx = 1
            for row in 1:(length(lib))
                for col in 1:row
                    values[idx] = obj_factor * Q[row, col]
                    idx += 1
                end
            end
        end
    end
                    
    return eval_f, eval_grad_f, eval_g, eval_jac_g, eval_h
end

function fit(lib::BinnedLibrary, data::KinematicData)
    nbins = length(data)
    n = length(lib) # number of variables
    x_L = zeros(n)
    x_U = ones(n)
    m = nbins + 1 # number of constraints
    g_L = vcat(zeros(nbins), 1)
    g_U = ones(m)
    nele_jac = m * n
    nele_hess = sum(1:n)
    eval_f, eval_grad_f, eval_g, eval_jac_g, eval_h = objective(lib, data)

    println("Creating problem")
    prob = createProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                         eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)
    prob.x = ones(n) / n
    println("Solving problem")
    status = solveProblem(prob)
    println(Ipopt.ApplicationReturnStatus[status])
    println("Weights: ", prob.x)
    println("Objective: ", prob.obj_val)
    return prob
end

end
