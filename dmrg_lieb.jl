using Plots
using ITensors
using Base.Threads
using LinearAlgebra
Threads.nthreads()
include("my_lattice.jl")
# 关闭BLAS和Strided多线程
BLAS.set_num_threads(1)
ITensors.Strided.disable_threads()
# 根据需要开启块稀疏多线程
ITensors.enable_threaded_blocksparse(true)
function main(; Nx=6, Ny=4, U=4.0, t=1.0, tp=0)
    N = 2Nx * 2Ny
    nsweeps = 10
    maxdim = [100, 200, 400, 800, 1600]
    cutoff = [1E-6]
    noise = [1E-6, 1E-7, 1E-8, 0.0]

    sites = siteinds("Electron", N;  conserve_qns=true)

    lattice1 = lieb_lattice_intra(Nx, Ny; yperiodic=true)
    lattice2 = lieb_lattice_inter(Nx, Ny; yperiodic=true)
    lattice3 = lieb_lattice_3NN(Nx, Ny; yperiodic=true)
    lattice4 = lieb_lattice(Nx, Ny; yperiodic=true)
    os = OpSum()
    for b in lattice1
        os += -t, "Cdagup", b.s1, "Cup", b.s2
        os += -t, "Cdagup", b.s2, "Cup", b.s1
        os += -t, "Cdagdn", b.s1, "Cdn", b.s2
        os += -t, "Cdagdn", b.s2, "Cdn", b.s1
    end
    for b in lattice2
        os += -t, "Cdagup", b.s1, "Cup", b.s2
        os += -t, "Cdagup", b.s2, "Cup", b.s1
        os += -t, "Cdagdn", b.s1, "Cdn", b.s2
        os += -t, "Cdagdn", b.s2, "Cdn", b.s1
    end
    for b in lattice3
        os += -tp, "Cdagup", b.s1, "Cup", b.s2
        os += -tp, "Cdagup", b.s2, "Cup", b.s1
        os += -tp, "Cdagdn", b.s1, "Cdn", b.s2
        os += -tp, "Cdagdn", b.s2, "Cdn", b.s1
    end
    for n in 1:N
        os += U, "Nupdn", n
    end
    H = MPO(os, sites)
    # Half filling
    # state = [(mod(n,4)==1 || mod(n,4) == 2) ? "Up" : "Dn" for n in 1:N]
    state = [isodd(n) ? "Dn" : "Up" for n in 1:N]

    # Initialize wavefunction to a random MPS
    # of bond-dimension 10 with same quantum
    # numbers as `state`
    psi0 = randomMPS(sites, state; linkdims=3)

    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise)
    @show t, U
    @show flux(psi)
    @show maxlinkdim(psi)
    @show energy
    return psi
end
# # 如果需要，在计算完成后恢复默认设置
# ITensors.enable_threaded_blocksparse(false)
psi = main()
magz = expect(psi, "Sz")
sum(magz)
plot_lieb_lattice_magnetic_moments(6,4,magz)
