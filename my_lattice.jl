using ITensors
function plot_lattice(lattice::Lattice)
  # Create the plot
  plt = plot(size=(600,600), axis_buffer=0.1, reuse=false)

  # Plot bonds as lines without a legend
  for bond in lattice
      plot!(plt, [bond.x1, bond.x2], [bond.y1, bond.y2], color=:black, label="")
  end

  # Plot sites as points
  scatter!(plt, [bond.x1 for bond in lattice], [bond.y1 for bond in lattice], 
           label="", color=:blue, markersize=5)
  scatter!(plt, [bond.x2 for bond in lattice], [bond.y2 for bond in lattice], 
          label="", color=:blue, markersize=5)

  # Configure the plot
  plot!(plt, xlabel="x", ylabel="y", title="Lieb Lattice")
  plot!(plt, axis=:equal) # Ensure the x and y axes have the same scale
  display(plt) # Display the plot
end
function lieb_lattice_intra(Nx::Int, Ny::Int; yperiodic=false)::Lattice
  yperiodic = yperiodic && (Ny > 2)
  N = 2Nx * 2Ny
  Nbond = 2 * Nx * Ny
  latt = Lattice(undef, Nbond)
  b = 0
  for n in 1:N
    x = div(n - 1, 2Ny) + 1
    y = mod(n - 1, 2Ny) + 1
    if mod(x, 2) == 1 && mod(y, 2) == 1
      latt[b += 1] = LatticeBond(n, n + 2Ny, x, y, x + 1, y)
      latt[b += 1] = LatticeBond(n, n + 1, x, y, x, y + 1)
    end
  end
  return latt
  
end

function lieb_lattice_inter(Nx::Int, Ny::Int; yperiodic=false)::Lattice
    yperiodic = yperiodic && (Ny > 2)
    N = 2Nx * 2Ny
    Nbond = 2 * Nx * Ny - Ny + (yperiodic ? 0 : -Nx)
    latt = Lattice(undef, Nbond)
    b = 0
    for n in 1:N
      x = div(n - 1, 2Ny) + 1
      y = mod(n - 1, 2Ny) + 1
      if x > 1 && mod(x, 2) == 1 && mod(y, 2) == 1
        latt[b += 1] = LatticeBond(n, n - 2Ny, x, y, x - 1, y)
      end
      if Ny > 1
        if y >1 && mod(x, 2) == 1 && mod(y, 2) == 1
          latt[b += 1] = LatticeBond(n, n - 1, x, y, x, y - 1)
        end
        if yperiodic && y == 1 && mod(x, 2) == 1
          latt[b += 1] = LatticeBond(n, n + 2Ny - 1, x, y, x, y + 2Ny - 1)
        end
      end
    end
    return latt
  end
  
  function lieb_lattice_3NN(Nx::Int, Ny::Int; yperiodic=false)::Lattice
    yperiodic = yperiodic && (Ny > 2)
    N = 2Nx * 2Ny
    Nbond = 2 * 3 * Nx * Ny - 3Ny + (yperiodic ? 0 : -3Nx)
    latt = Lattice(undef, Nbond)
    b = 0
    for n in 1:N
      x = div(n - 1, 2Ny) + 1
      y = mod(n - 1, 2Ny) + 1
      if x < 2Nx - 1 && (mod(y, 2) == 1 || mod(x, 2) == 1)
        latt[b += 1] = LatticeBond(n, n + 4Ny, x, y, x + 2, y)
      end
      if Ny > 1
        if y < 2Ny - 1 && (mod(y, 2) == 1 || mod(x, 2) == 1)
          latt[b += 1] = LatticeBond(n, n + 2, x, y, x, y + 2)
        end
        if yperiodic && y == 2Ny - 1 && (mod(y, 2) == 1 || mod(x, 2) == 1)
          latt[b += 1] = LatticeBond(n, n - 2Ny + 2, x, y, x, y - 2Ny + 2)
        end
        if yperiodic && y == 2Ny && (mod(y, 2) == 1 || mod(x, 2) == 1)
          latt[b += 1] = LatticeBond(n, n - 2Ny + 2, x, y, x, y - 2Ny + 2)
        end
      end
    end
    return latt
  end