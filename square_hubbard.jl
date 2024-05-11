using ITensors
using Plots
include("my_lattice.jl")

function main(; Nx=8, Ny=4, U=50, t=1.0)
  N = Nx * Ny

  nsweeps = 20
  maxdim = [100, 200, 400, 800, 1600]
  cutoff = [1E-6]
  noise = [1E-6, 1E-7, 1E-8, 0.0]

  sites = siteinds("Electron", N; conserve_nf=true)

  lattice = square_lattice(Nx, Ny; yperiodic=true)

  os = OpSum()
  for b in lattice
    os += -t, "Cdagup", b.s1, "Cup", b.s2
    os += -t, "Cdagup", b.s2, "Cup", b.s1
    os += -t, "Cdagdn", b.s1, "Cdn", b.s2
    os += -t, "Cdagdn", b.s2, "Cdn", b.s1
  end
  for n in 1:N
    os += U, "Nupdn", n
  end
  H = MPO(os, sites)

  # Half filling
  state = [isodd(n) ? "Dn" : "Up" for n in 1:N]

  # Initialize wavefunction to a random MPS
  # of bond-dimension 10 with same quantum
  # numbers as `state`
  psi0 = randomMPS(sites, state; linkdims=5)

  energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise)
  @show t, U
  @show flux(psi)
  @show maxlinkdim(psi)
  @show energy

  return psi
end

psi = main()

magz = expect(psi,"Sz")
sum(magz)
plot_square_lattice_magnetic_moments(4, 2, magz)
