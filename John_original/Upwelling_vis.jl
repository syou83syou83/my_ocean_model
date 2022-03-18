# ## Visualize the results

using JLD2, Plots
using Printf
using Oceananigans
using Oceananigans.Units: hours, hour, day, days, minute, minutes

# Re-create the model grid to get the coordinate arrays
Lx=200e3;
Ly=200e3;
Lz=2e3;
Nx=128;
Ny=128;
Nz=32;
# Vertical grid stretching parameters
refinement = 4 # controls spacing near surface (higher means finer spaced)
stretching = 10   # controls rate of stretching at bottom 
## Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz
## Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement
## Bottom-intensified stretching function 
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
## Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)
grid = RectilinearGrid(size = (Nx, Ny, Nz), 
                                         x = (0, Lx),
                                         y = (0, Ly),
                                         z = z_faces,
                                          topology = (Bounded, Periodic, Bounded))
model = NonhydrostaticModel(
        #   architecture = CPU(),
                   grid = grid,
            #   advection = UpwindBiasedFifthOrder(),
            # timestepper = :RungeKutta3,
            #    coriolis = coriolis,
                tracers = (:T, :S, :C),
            #    buoyancy = buoyancy,
            #    particles = lagrangian_particles,
                # closure = (Laplacian_vertical_diffusivity, biharmonic_horizontal_diffusivity),
                # forcing = (C=C_forcing,u=IB_u_forcing, v=IB_v_forcing, w=IB_w_forcing),
    # boundary_conditions = (u=u_bcs, v=v_bcs)
)

## Coordinate arrays
xw, yw, zw = nodes(model.velocities.w)
xu, yu, zu = nodes(model.velocities.u)
xv, yv, zv = nodes(model.velocities.v)

xT, yT, zT = nodes(model.tracers.T)
xS, yS, zS = nodes(model.tracers.S)
xC, yC, zC = nodes(model.tracers.C)

## Open the file with our data
file_xz = jldopen("./slice_xz.jld2")
file_xy = jldopen("./slice_xy.jld2")
# file_particles = jldopen(simulation.output_writers[:particles].filepath)

## Extract a vector of iterations
iterations = parse.(Int, keys(file_xy["timeseries/t"]))

""" Returns colorbar levels equispaced between `(-clim, clim)` and encompassing the extrema of `c`. """
function divergent_levels(c, clim, nlevels=21)
    cmax = maximum(abs, c)
    levels = clim > cmax ? range(-clim, stop=clim, length=nlevels) : range(-cmax, stop=cmax, length=nlevels)
    return (levels[1], levels[end]), levels
end

""" Returns colorbar levels equispaced between `clims` and encompassing the extrema of `c`."""
function sequential_levels(c, clims, nlevels=20)
    levels = range(clims[1], stop=clims[2], length=nlevels)
    cmin, cmax = minimum(c), maximum(c)
    cmin < clims[1] && (levels = vcat([cmin], levels))
    cmax > clims[2] && (levels = vcat(levels, [cmax]))
    return clims, levels
end
nothing # hide

times = [file_xy["timeseries/t/$iter"] for iter in iterations]
intro = searchsortedfirst(times, 0minutes)

anim = @animate for (i, iter) in enumerate(iterations[intro:end])

    @info "Drawing frame $i from iteration $iter..."

    t = file_xy["timeseries/t/$iter"]
    u_xy = file_xy["timeseries/u/$iter"][:, :, 1]
    v_xy = file_xy["timeseries/v/$iter"][:, :, 1]
    w_xy = file_xy["timeseries/w/$iter"][:, :, 1]
    u_xz = file_xz["timeseries/u/$iter"][:, 1, :]
    v_xz = file_xz["timeseries/v/$iter"][:, 1, :]
    w_xz = file_xz["timeseries/w/$iter"][:, 1, :]

    T_xz = file_xz["timeseries/T/$iter"][:, 1, :]
    S_xz = file_xz["timeseries/S/$iter"][:, 1, :]
    C_xz = file_xz["timeseries/C/$iter"][:, 1, :]
    T_xy = file_xy["timeseries/T/$iter"][:, :, 1]
    S_xy = file_xy["timeseries/S/$iter"][:, :, 1]
    C_xy = file_xy["timeseries/C/$iter"][:, :, 1]

    # particles = file_particles["timeseries/particles/$iter"]

    ulims, ulevels = divergent_levels(u_xz,1)
    wlims, wlevels = divergent_levels(w_xz, 2e-2)
    Tlims, Tlevels = sequential_levels(T_xy, (15, 20))
    Slims, Slevels = sequential_levels(S_xz, (33.75, 35.5))
    Clims, Clevels = sequential_levels(C_xy, (0, 1))

    kwargs_xz = (linewidth=0, xlabel="x (m)", ylabel="z (m)", aspectratio=grid.Lx/grid.Lz,
              xlims=(0, grid.Lx), ylims=(-grid.Lz, 0))
    kwargs_xy = (linewidth=0, xlabel="x (km)", ylabel="y (km)", aspectratio=1,
              xlims=(0, grid.Lx/1e3), ylims=(0, grid.Ly/1e3))

    w_plot = contourf(xw, yw, w_xy'; color=:balance, clims=wlims, levels=wlevels, kwargs_xy...)
    T_xz_plot = contourf(xT, zT, T_xz'; color=:thermal, clims=Tlims, levels=Tlevels, kwargs_xz...)
    # scatter!(particles.x,particles.z, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    T_xy_plot = contourf(xT/1e3, yT/1e3, T_xy'; color=:thermal, clims=Tlims, levels=Tlevels, kwargs_xy...)
    # scatter!(particles.x,particles.y, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    S_xz_plot = contourf(xS, zS, S_xz'; color=:thermal, clims=Slims, levels=Slevels, kwargs_xz...)
    # scatter!(particles.x,particles.z, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    S_xy_plot = contourf(xS/1e3, yS/1e3, S_xy'; color=:thermal, clims=Slims, levels=Slevels, kwargs_xy...)
    # scatter!(particles.x,particles.y, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    C_xz_plot = heatmap(xC, zC, C_xz'; color=:thermal, clims=(0,10), levels=Clevels, kwargs_xz...)
    # scatter!(particles.x,particles.z, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    C_xy_plot = contourf(xC/1e3, yC/1e3, C_xy'; color=:haline, clims=(0,1), levels=Clevels, kwargs_xy...)
    # scatter!(particles.x,particles.y, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    u_xy_plot = heatmap(xu, yu, u_xy'; color=:balance, clims=ulims, levels=ulevels, kwargs_xy...)    
    v_xy_plot = heatmap(xv, yv, v_xy'; color=:balance, clims=ulims, levels=ulevels, kwargs_xy...)    
    u_xz_plot = heatmap(xu, zu, u_xz'; color=:balance, clims=ulims, levels=ulevels, kwargs_xz...)    
    v_xz_plot = heatmap(xv, zv, v_xz'; color=:balance, clims=ulims, levels=ulevels, kwargs_xz...)   

    w_title = @sprintf("vertical velocity (m s⁻¹), t = %s", prettytime(t))
    T_title = "temperature (ᵒC)"
    S_title = "salinity (g kg⁻¹)"
    u_title = "u (m/s)"
    v_title = "v (m/s)"
    C_title = "tracer concentration"
    #ν_title = "eddy viscosity (m² s⁻¹)"

    ## Arrange the plots side-by-side.
#    plot(T_xy_plot, C_xy_plot, u_xz_plot, C_xz_plot, layout=(2, 2), size=(1200, 1200),
#         title=[T_title C_title T_title C_title])
    plot(T_xy_plot, C_xy_plot, size=(1200,600), title=[T_title C_title])

    #iter == iterations[end] && close(file)
end

mp4(anim, "Upwelling.mp4", fps = 20) # hide
