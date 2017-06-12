using Schwarz: Constants, PowerLawPotential, OrbitLibrary
c = Constants

r0 = 10 * c.kpc
rho0 = 4e7 * c.Msun / c.kpc^3
alpha = -2.5
phi = PowerLawPotential(r0, rho0, alpha)

Nmax = 16
reps = 5
Nv = 4
rmax = 40 * c.kpc
multithreaded = true
n = nworkers()

sizes = []
means = []
stds = []

println("Pre-compiling OrbitLibrary")
lib = OrbitLibrary(1, 1, rmax, phi; multithreaded = true)
Nr = 1
lib, elapsed = @timed OrbitLibrary(Nr, Nv, rmax, phi; multithreaded = multithreaded)

for Nr in 1:Nmax
    println("Timing OrbitLibrary")
    println("    Nr = $Nr, Nv = $Nv")
    times = Array{Float64, 1}(reps)
    for i in 1:reps
        lib, elapsed = @timed OrbitLibrary(Nr, Nv, rmax, phi; multithreaded = multithreaded)
        times[i] = elapsed
    end
    push!(sizes, Nr)
    push!(means, mean(times))
    push!(stds, std(times))
end

writedlm("orbitLibrary_$n.txt", [sizes means stds])
