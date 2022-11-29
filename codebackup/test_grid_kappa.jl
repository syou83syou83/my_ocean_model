
using Oceananigans
using Plots
Lx = 1  #20
Ly = 1  #20
Nx = 1
Ny = 1
Nz = 50 # number of points in the vertical direction
Lz = 600 # domain depth

# Generate vertically stretched grid 
refinement = 10 # controls spacing near surface (higher means finer spaced)  #10         1.2
stretching = 0.5   # controls rate of stretching at bottom                 #5.754           5
# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz
# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement
# Bottom-intensified stretching function 
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)
grid = RectilinearGrid(size = (Nx, Ny, Nz), 
                       x = (0, Lx),
                       y = (0, Ly),
                       z = z_faces)        #z = z_faces

Plots.plot(grid.Δzᵃᵃᶜ[1:grid.Nz], grid.zᵃᵃᶜ[1:grid.Nz],
                       marker = :circle,
                       ylabel = "Depth (m)",
                       xlabel = "Vertical spacing (m)",
                       legend = nothing)





k(z) = 8e-2*max(1-(z+100/2)^2/(100/2)^2,0)+1e-4; #setup viscosity and diffusivity in the following Model instantiation
k2(z) = 8e-2*max(1-(z+100/2)^2/(100/2)^2,0)+1e-4; #setup viscosity and diffusivity in the following Model instantiation

w=-600:0
Plots.plot(k2.(w),w)

mld=100
k_mld=8e-2
k_interior=1e-4
a=(k_mld-k_interior)/mld^2
b=2*a*mld

# k3(z)=a*z^2+b*z+k_mld
k3(z)=a*max(z,-mld)^2+b*max(z,-mld)+k_mld

Plots.plot(k3.(w),w)
