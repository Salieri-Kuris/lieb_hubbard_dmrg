using ITensors
function plot_bond(lattice::Lattice)
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

function lieb_lattice(Nx::Int, Ny::Int; yperiodic=false)::Lattice
  yperiodic = yperiodic && (Ny > 1)
  N = 2Nx * 2Ny
  Nbond = 2 * Nx * Ny + 2 * Nx * Ny - Ny + (yperiodic ? 0 : -Nx)
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
# function lieb_lattice_intra(Nx::Int, Ny::Int; yperiodic=false)::Lattice
#   yperiodic = yperiodic && (Ny > 2)
#   N = 2Nx * 2Ny
#   Nbond = 2 * Nx * Ny
#   latt = Lattice(undef, Nbond)
#   b = 0
#   for n in 1:N
#     x = div(n - 1, 2Ny) + 1
#     y = mod(n - 1, 2Ny) + 1
#     if mod(x, 2) == 1 && mod(y, 2) == 1
#       latt[b += 1] = LatticeBond(n, n + 2Ny, x, y, x + 1, y)
#       latt[b += 1] = LatticeBond(n, n + 1, x, y, x, y + 1)
#     end
#   end
#     return latt
  
# end
function lieb_lattice_intra(Nx::Int, Ny::Int; yperiodic=false)::Lattice
  yperiodic = yperiodic && (Ny > 2)
  N = 3 * Nx * Ny
  Nbond = 2 * Nx * Ny
  latt = Lattice(undef, Nbond)
  b = 0
  for n in 1:N
    x = 2 * div(n, 3Ny) + ((mod(n , 3Ny) != 0) ? (mod(n, 3Ny) > 2Ny ? 2 : 1) : 0)
    y = (mod(n, 3Ny) > 2Ny  ? 2(mod(n, 3Ny) - 2Ny)-1 : (mod(n, 3Ny) == 0 ? 2Ny-1 : mod(n, 3Ny)))
    if mod(x, 2) == 1 && mod(y, 2) == 1
      latt[b += 1] = LatticeBond(n, n + 2Ny, x, y, x + 1, y)
      latt[b += 1] = LatticeBond(n, n + 1, x, y, x, y + 1)
    end
  end
    return latt
  
end

# function lieb_lattice_inter(Nx::Int, Ny::Int; yperiodic=false)::Lattice
#     yperiodic = yperiodic && (Ny > 2)
#     N = 2Nx * 2Ny
#     Nbond = 2 * Nx * Ny - Ny + (yperiodic ? 0 : -Nx)
#     latt = Lattice(undef, Nbond)
#     b = 0
#     for n in 1:N
#       x = div(n - 1, 2Ny) + 1
#       y = mod(n - 1, 2Ny) + 1
#       if x > 1 && mod(x, 2) == 1 && mod(y, 2) == 1
#         latt[b += 1] = LatticeBond(n, n - 2Ny, x, y, x - 1, y)
#       end
#       if Ny > 1
#         if y >1 && mod(x, 2) == 1 && mod(y, 2) == 1
#           latt[b += 1] = LatticeBond(n, n - 1, x, y, x, y - 1)
#         end
#         if yperiodic && y == 1 && mod(x, 2) == 1
#           latt[b += 1] = LatticeBond(n, n + 2Ny - 1, x, y, x, y + 2Ny - 1)
#         end
#       end
#     end
#     return latt
#   end

function lieb_lattice_inter(Nx::Int, Ny::Int; yperiodic=false)::Lattice
  yperiodic = yperiodic && (Ny > 2)
  N = 3 * Nx * Ny
  Nbond = 2 * Nx * Ny - Ny + (yperiodic ? 0 : -Nx)
  latt = Lattice(undef, Nbond)
  b = 0
  for n in 1:N
    x = 2 * div(n, 3Ny) + ((mod(n , 3Ny) != 0) ? (mod(n, 3Ny) > 2Ny ? 2 : 1) : 0)
    y = (mod(n, 3Ny) > 2Ny  ? 2(mod(n, 3Ny) - 2Ny)-1 : (mod(n, 3Ny) == 0 ? 2Ny-1 : mod(n, 3Ny)))
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

function plot_lieb_lattice_magnetic_moments(Nx, Ny, magz)
  @assert length(magz) == 4 * Nx * Ny "The length of magnetic moments vector must equal 4 * Nx * Ny."
  plot_size = max(Nx, Ny) * 100
  # Create a scatter plot for the magnetic moments
  plt = scatter(
      # x coordinates for each site
      [div(n - 1, 2 * Ny) + 1 for n in 1:length(magz)],
      # y coordinates for each site
      [mod(n - 1, 2 * Ny) + 1 for n in 1:length(magz)],
      # Magnetic moment values as annotations
      series_annotations=[(string(round(val, digits=4)), font(8,colorant"red")) for val in magz],
      # Markersize is set to a small value to make sure the annotations are visible
      markersize=6,
      size=(plot_size, plot_size),
      # Set the aspect ratio to 1 for a square lattice
      aspect_ratio=0.6,
      # Optionally, set a title for the plot
      title="Magnetic Moments on Square Lattice",
      # Set the x and y limits to match the lattice size
      xlims=(1, 2 * Nx),
      ylims=(1, 2 * Ny),
      # Other styling options
      dpi=1000,
      label = ""
  )

  # Set the xticks and yticks to match the coordinates
  xticks!(plt, 1:2:Nx*2-1)
  yticks!(plt, 1:2:Ny*2-1)
  savefig(plt, "lieb_mag")
  # Display the plot
  display(plt)
  return plt
end

