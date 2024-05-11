using Plots
using ITensors
using Base.Threads
using LinearAlgebra
using Formatting
Threads.nthreads()
include("my_lattice.jl")
# 关闭BLAS和Strided多线程
BLAS.set_num_threads(1)
ITensors.Strided.disable_threads()
# 根据需要开启块稀疏多线程
ITensors.enable_threaded_blocksparse(true)
Nx = 8
Ny = 5
U = 400.0
del = 0.0
tp = 0.0

function main(; Nx, Ny, U, del, tp)
    N = 3 * Nx * Ny
    t = 1.0
    nsweeps = 20
    maxdim = [100, 200, 400, 800, 1600]
    cutoff = [1E-6]
    noise = [1E-6, 1E-7, 1E-8, 0.0]

    sites = siteinds("Electron", N;  conserve_nf=true)

    lattice1 = lieb_lattice_intra(Nx, Ny; yperiodic=true)
    lattice2 = lieb_lattice_inter(Nx, Ny; yperiodic=true)
    lattice3 = lieb_lattice_3NN(Nx, Ny; yperiodic=true)
    os = OpSum()
    t1 =  t * (1 + del)
    t2 =  t * (1 - del)
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
    # for b in lattice3
    #     os += -tp, "Cdagup", b.s1, "Cup", b.s2
    #     os += -tp, "Cdagup", b.s2, "Cup", b.s1
    #     os += -tp, "Cdagdn", b.s1, "Cdn", b.s2
    #     os += -tp, "Cdagdn", b.s2, "Cdn", b.s1
    # end
    for n in 1:N
        os += U, "Nupdn", n
    end
    H = MPO(os, sites)
    # Half filling
    # state = [(mod(n,4)==1 || mod(n,4) == 2) ? "Up" : "Dn" for n in 1:N]
    # state = [mod(n, 3)>1 ? "Dn" : "Up" for n in 1:N]
    state = [(mod(2 * div(n, 3Ny) + ((mod(n , 3Ny) != 0) ? (mod(n, 3Ny) > 2Ny ? 2 : 1) : 0), 2) == 1 && mod((mod(n, 3Ny) > 2Ny  ? 2(mod(n, 3Ny) - 2Ny)-1 : (mod(n, 3Ny) == 0 ? 2Ny-1 : mod(n, 3Ny))), 2) == 1) ? "Dn" : "Up" for n in 1:N]

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
# # 如果需要，在计算完成后恢复默认设置
# ITensors.enable_threaded_blocksparse(false)
psi = main(;Nx, Ny, U, del, tp)
magz = expect(psi, "Sz")
sum(magz)
plot_lieb_lattice_magnetic_moments(Nx, Ny, U, del, tp, magz)
anti_order(magz, 4)
lattice1 = lieb_lattice_intra(Nx, Ny; yperiodic=true)
lattice2 = lieb_lattice_inter(Nx, Ny; )
plot_bond(lattice1)
t = 1.0
os = OpSum()
    t1 =  t * (1 + del)
    t2 =  t * (1 - del)
    for b in lattice1
        os += -t1, "Cdagup", b.s1, "Cup", b.s2
        os += -t1, "Cdagup", b.s2, "Cup", b.s1
        os += -t1, "Cdagdn", b.s1, "Cdn", b.s2
        os += -t1, "Cdagdn", b.s2, "Cdn", b.s1
    end
    for b in lattice2
        os += -t2, "Cdagup", b.s1, "Cup", b.s2
        os += -t2, "Cdagup", b.s2, "Cup", b.s1
        os += -t2, "Cdagdn", b.s1, "Cdn", b.s2
        os += -t2, "Cdagdn", b.s2, "Cdn", b.s1
    end
os
sites = siteinds("Electron", 72; conserve_qns=true)
state = [mod(n, 3)>1 ? "Dn" : "Up" for n in 1:72]
psi0 = randomMPS(sites, state; linkdims=5)
mag_random = expect(psi0, "Sz")
sum(mag_random)