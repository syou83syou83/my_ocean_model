"""
Nonhydrostaticmodel, with buoyancy as tracer(:b) and particles(:c), save all data except halo regions 
10.3.2022
Si Chen

"""

using Random
using Printf
using Plots
using JLD2
using Statistics

using Oceananigans
using Oceananigans.Advection
using Oceananigans.Operators
using Oceananigans.Fields: ZeroField
using Oceananigans.Fields: BackgroundFields, Field, tracernames, VelocityFields, TracerFields, PressureFields
using Oceananigans.Units: minute, minutes, hour, hours

Nz = 24          # number of points in the vertical direction
Lz = 32          # (m) domain depth
Nx = 32
Ny = 32

refinement = 1.2 # controls spacing near surface (higher means finer spaced)
stretching = 12  # controls rate of stretching at bottom

# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz

# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement

# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(size = (Nx, Ny, Nz),
                          x = (0, Nx*2),
                          y = (0, Ny*2),
                          z = z_faces)

plot(grid.Δzᵃᵃᶜ[1:grid.Nz], grid.zᵃᵃᶜ[1:grid.Nz],
     marker = :circle,
     ylabel = "Depth (m)",
     xlabel = "Vertical spacing (m)",
     legend = nothing)

grid.zᵃᵃᶜ;



u₁₀ = 10    # m s⁻¹, average wind velocity 10 meters above the ocean
cᴰ = 2.5e-3 # dimensionless drag coefficient
ρₒ = 1026 # kg m⁻³, average density at the surface of the world ocean
ρₐ = 1.225  # kg m⁻³, average density of air at sea-level

Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²   kinematic stress: stress divided by density 

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))

B₀ = -4.24e-8    #m²s⁻³
N² = 9e-6    #dbdz=N^2, s⁻²
buoyancy_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(B₀),
                                       bottom = GradientBoundaryCondition(N²))

# Add an extra advection term to the tracer equation to account for a 'slip' velocity

ws = 0.001;   # slip velocity (in m/s)
λ = 2;      # 2 decay length of slip velocity at z=0 (m)，Si:deeper than this, slip velocity is very close to ws.

# create an empty array to store the slip velocity (this is a tuple)
slip_vel = (u=zeros(0:Nx+2,0:Ny+2,0:Nz+2),v=zeros(0:Nx+2,0:Ny+2,0:Nz+2),w=zeros(0:Nx+2,0:Ny+2,0:Nz+2))
# because div_Uc is 3D, so slip_vel has to be defined 3D.
# zeros function in Oceananigans is different in Julia because in Julia zeros(0:2) erros.
# to learn the dimensions and indices see the following  
# when build it has with 1×1×1 halo, after calculation, it has 3×3×3 halo, actuall 3 in both sides of each direction  

# define the non-zero components(s) of the slip velocity
# to conserve mass, the slip velocity should =0 at z=0
# The width of the tanh will depend on the advection scheme and integrated tracer conservation should be checked
for k=0:Nz+2
  slip_vel.w[:,:,k].+=ws*tanh(max(-grid.zᵃᵃᶠ[k]/λ,0));   # don't need .+=  this "+"
end

# define the RHS forcing term for advection by the slip velocity using the inbuilt advection operator
# Different advection schemes can be specified here

#slip_advection(i, j, k, grid, clock, model_fields) = - div_Uc(i,j,k,grid,CenteredSecondOrder(),slip_vel,model_fields.c)  
slip_advection(i, j, k, grid, clock, model_fields) = - div_Uc(i,j,k,grid,WENO5(grid = grid, stretched_smoothness = true),slip_vel,model_fields.c)  



# because div_Uc is 3D, so slip_vel has to be defined 3D. 
# see: /src/Advection/centered_second_order.jl here it also defined a constant: const C2 = CenteredSecondOrder L14
# try different scheme like WENO5 https://github.com/CliMA/Oceananigans.jl/blob/main/src/Advection/weno_fifth_order.jl

slip_forcing = Forcing(slip_advection, discrete_form=true)

#see: https://clima.github.io/OceananigansDocumentation/stable/model_setup/forcing_functions/
    # 1 Forcing functions with external parameters means only influencing x,y,z,t
    # 2 Forcing functions that depend on model fields means may influence u,v,w,S,T, etc.  
    # 3 "Discrete form" forcing functions means by using model_fields argument, it can access all the model field data u,v,w, 
         # and all tracers through like model_fields.b, model_fields.u, model_fields.c  , etc. 
    # "Discrete form" forcing functions




plot(slip_vel.w[1,1,1:Nz+1],grid.zᵃᵃᶠ[1:Nz+1],     
    marker = :circle,
    ylabel = "Depth (m)",
    xlabel = "Slip Velocities (ms⁻¹)",)

grid.zᵃᵃᶠ

model = NonhydrostaticModel(
                            advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            grid = grid,
                            tracers = (:b,:c),
                            coriolis = FPlane(f=1e-4),
                            buoyancy = BuoyancyTracer(),
                            closure = AnisotropicMinimumDissipation(),
                            forcing = (c=slip_forcing,),
                            boundary_conditions = (u=u_bcs, b=buoyancy_bcs))

model.forcing.c

# Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

# Velocity initial condition: random noise scaled by the friction velocity.
uᵢ(x, y, z) = sqrt(abs(Qᵘ)) * 1e-3 * Ξ(z)

cᵢ(x, y, z) = exp(z/10)



#set!(model, u=uᵢ, w=uᵢ, b=1e-6)

mixed_layer_depth = -20  # m
∂b∂z(z) = z > mixed_layer_depth ? 0 : 9e-6   # https://juliacn.gitlab.io/JuliaZH.jl/manual/functions.html

b₀(x, y, z) = ∂b∂z(z) * z 
#https://clima.github.io/OceananigansDocumentation/stable/generated/convecting_plankton/

set!(model, u=uᵢ, w=uᵢ, b=b₀, c=cᵢ)


simulation = Simulation(model, Δt=10.0, stop_iteration=2)  # Simulation(model, Δt=10.0, stop_time=60minutes)

wizard = TimeStepWizard(cfl=1.0, max_change=1.1, max_Δt=1minute)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))


"""
julia> simulation = Simulation(model, Δt=0.9, stop_iteration=100)
julia> wizard = TimeStepWizard(cfl=0.2)
julia> simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(4))

Then when `run!(simulation)` is invoked, the time-step `simulation.Δt` will be updated every 4 iterations.
Note that the name `:wizard` is unimportant.
"""


# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim),
#src/Simulations/simulation.jl#L111-L115
                                prettytime(sim),
                                prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.w),
                                prettytime(sim.run_wall_time))
#prettytime: Convert a floating point value to a more human-friendly formatted string with three decimal places.
# also a generic function 
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(10))

# Create a NamedTuple with eddy viscosity
eddy_viscosity = (; νₑ = model.diffusivity_fields.νₑ)

simulation.output_writers[:slices] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_viscosity),
                           prefix = "ocean_wind_mixing_and_convection_BuoyancyPlusParticleW5"*string(λ),
                     #field_slicer = FieldSlicer(j=Int(grid.Ny/2)), #here save only along y direction
                        field_slicer = FieldSlicer(),   #include all data except halo retions
                          # indices = (:, 16, :), #error 
    
# field_slicer = nothing` means no slicing occurs, so that all field data, including halo regions, is saved.
#Default: FieldSlicer(), which slices halo regions. 
#FieldSlicer() Si: only grid includes halos,+-3,and u,v,w, and all field data don't include halo. 
#don't forget z include lower bound but face and center have different dimensions.
                         schedule = TimeInterval(1minute),
                            force = true)


run!(simulation)   # Δt is changing.

""" Returns colorbar levels equispaced between `(-clim, clim)` and encompassing the extrema of `c`. """
function divergent_levels(c, clim, nlevels=21) 
    # the input clim is only one number,divergent_levels is for those the output is 
    # symmetric(on positive and negative sides), like velocity
    cmax = maximum(abs, c)  
    #Compute the maximum value by calling the function abs on each element of c over the given dimensions.
    # i.e. calculate abs first than max
    levels = clim > cmax ? range(-clim, stop=clim, length=nlevels) : range(-cmax, stop=cmax, length=nlevels)
    return (levels[1], levels[end]), levels
end

""" Returns colorbar levels equispaced between `clims` and encompassing the extrema of `c`."""
function sequential_levels(c, clims, nlevels=20)
        # the input clims is a range, sequential_levels is for those the output is 
        # nonsymmetric(on positive and negative sides), like scalars 
    levels = range(clims[1], stop=clims[2], length=nlevels)
    cmin, cmax = minimum(c), maximum(c)
    cmin < clims[1] && (levels = vcat([cmin], levels)) 
    cmax > clims[2] && (levels = vcat(levels, [cmax]))
    return (levels[1], levels[end]), levels    # I think the return should like this
    #return clims, levels
end

# Coordinate arrays
xw, yw, zw = nodes(model.velocities.w)  #there is no real practical data in “model”, just the structure
xb, yb, zb = nodes(model.tracers.b)   # velocity field is located on Face while scalar is located on Center

# Open the file with our data
file = jldopen(simulation.output_writers[:slices].filepath)

# Extract a vector of iterations
iterations = parse.(Int, keys(file["timeseries/t"])) # turn an array of strings into numbers 

times = [file["timeseries/t/$iter"] for iter in iterations]
intro = searchsortedfirst(times, 10minutes);  
# show results from 10 minutes 
# https://cn.julialang.org/JuliaZH.jl/latest/base/sort/#Base.Sort.searchsortedfirst
keys(file)

anim = @animate for (i, iter) in enumerate(iterations[intro:end])

    @info "Drawing frame $i from iteration $iter..."

    t = file["timeseries/t/$iter"]
    w = file["timeseries/w/$iter"][:, 1, :]
    b = file["timeseries/b/$iter"][:, 1, :]
    c = file["timeseries/c/$iter"][:, 1, :]
    #S = file["timeseries/S/$iter"][:, 1, :]
    νₑ = file["timeseries/νₑ/$iter"][:, 1, :]

    wlims, wlevels = divergent_levels(w, 2e-2)              #functions defined above
    blims, blevels = sequential_levels(b, (-3e-4, 5e-5))   # (-3e-4, 5e-5)
    clims, clevels = sequential_levels(c, (0, 1))
    #Slims, Slevels = sequential_levels(S, (35, 35.005))
    νlims, νlevels = sequential_levels(νₑ, (1e-6, 5e-3))

    kwargs = (linewidth=0, xlabel="x (m)", ylabel="z (m)", aspectratio=1,
              xlims=(0, grid.Lx), ylims=(-grid.Lz, 0))

    w_plot = contourf(xw, zw, w'; color=:balance, clims=wlims, levels=wlevels, kwargs...)
    b_plot = contourf(xb, zb, b'; color=:thermal, clims=blims, levels=blevels, kwargs...)
    c_plot = contourf(xb, zb, c'; color=:thermal, clims=clims, levels=clevels, kwargs...)
    # argument clims:  Fixes the limits of the colorbar.
        #https://docs.juliaplots.org/latest/generated/attributes_subplot/
    # argument levels: Determines contour levels for a contour type.
        #https://docs.juliaplots.org/latest/generated/attributes_series/ 
    
    #S_plot = contourf(xT, zT, S'; color=:haline,  clims=Slims, levels=Slevels, kwargs...)

    # We use a heatmap for the eddy viscosity to observe how it varies on the grid scale.
    ν_plot = heatmap(xb, zb, νₑ'; color=:thermal, clims=νlims, levels=νlevels, kwargs...)

    w_title = @sprintf("vertical velocity (ms⁻¹), t = %s", prettytime(t))
    b_title = "buoyancy (ms⁻²)"
    c_title = "tracer"
    #S_title = "salinity (g kg⁻¹)"
    ν_title = "eddy viscosity (m²s⁻¹)"

    # Arrange the plots side-by-side.
    plot(w_plot, b_plot, c_plot, ν_plot, layout=(2, 2), size=(1200, 600),
         title=[w_title b_title c_title ν_title])

    #iter == iterations[end] && close(file) # keep file open until all process is done
end

gif(anim, "ocean_wind_mixing_and_convection_BuoyancyPlusParticle.mp4", fps = 5)  # https://docs.juliaplots.org/latest/animations/

#file = jldopen(simulation.output_writers[:slices].filepath)
#iterations = parse.(Int, keys(file["timeseries/t"]));
#times = [file["timeseries/t/$iter"] for iter in iterations]
#intro = searchsortedfirst(times, 10minutes) # show results from 10 minutes 

# https://cn.julialang.org/JuliaZH.jl/latest/base/sort/#Base.Sort.searchsortedfirst
v_cell=zeros(1,Nz)          #cell volume 

for i in 1:Nz
    v_cell[i]=grid.Δxᶜᵃᵃ * grid.Δyᵃᶜᵃ * grid.Δzᵃᵃᶜ[i]  # calculate cell volume 1*Nz dimension
end


t=zeros(1,length(iterations)-intro+1)
c=zeros(1,length(iterations)-intro+1)

for (i, iter) in enumerate(iterations[intro:end])
    #println((i,iter))
    t[i] = file["timeseries/t/$iter"]        
    concentration=zeros(Ny,Nz)
    for j in 1:Nx
        concentration+=file["timeseries/c/$iter"][j,:,:].*v_cell
        #add volume-integrated concentration in x plane first
    end
    c[i] = sum(concentration)  # to obtain volume-integrated concentration in whole 3D domain in one timestep

end


plot(t[:]/60,c[:],ylims=[30000,40000],   
     marker = :circle,
     ylabel = "Volume-integrated concentration",
     xlabel = "Time (minutes)",
     #label = @printf("Max: %.5d, Min: %d, δ: %.6e ", maximum(c),minimum(c),(maximum(c)-minimum(c))/minimum(c) ))
     label = " Max=$(maximum(c)) \n Min=$(minimum(c)) \n δ=$((maximum(c)-minimum(c))/minimum(c)) ")

v_cell

file["timeseries/c/0"][1,:,:]

file["timeseries/c/0"][1,:,:].*v_cell

grid.Δzᵃᵃᶜ

model.velocities.u;

anim2 = @animate for (i, iter) in enumerate(iterations[intro:end])

    #@info "Drawing frame $i from iteration $iter..."

    t = file["timeseries/t/$iter"]
    u_end=file["timeseries/u/$iter"]
    v_end=file["timeseries/v/$iter"]
    w_end=file["timeseries/w/$iter"]
    
    u_avg=mean(u_end,dims=(1,2))
    v_avg=mean(v_end,dims=(1,2))
    w_avg=mean(w_end,dims=(1,2))
    
    u_prime=zeros(Nx,Ny,Nz)
    v_prime=zeros(Nx,Ny,Nz)
    w_prime=zeros(Nx,Ny,Nz)  #its dimension is actually Nx Ny Nz+1
    
    for k in 1:Nz
        u_prime[:,:,k]=u_end[:,:,k].-u_avg[k]  #broadcast array.-scalar 
        v_prime[:,:,k]=v_end[:,:,k].-v_avg[k]
        w_prime[:,:,k]=w_end[:,:,k].-w_avg[k]
    end
    u_rms=(mean((u_end-u_prime).^2,dims=(1,2))).^0.5
    v_rms=(mean((v_end-v_prime).^2,dims=(1,2))).^0.5
    w_rms=(mean((w_end[:,:,1:Nz]-v_prime).^2,dims=(1,2))).^0.5;
    
    plot(u_rms[:],model.grid.zᵃᵃᶠ[1:Nz],label = "uᵣₘₛ",
        ylabel = "Depth (m)",
        xlabel = "Velocities (ms⁻¹)",
        xlims = (0,0.3),
        title = @sprintf("RMS velocities (ms⁻¹), t = %s", prettytime(t)))    
    plot!(v_rms[:],model.grid.zᵃᵃᶠ[1:Nz],label = "vᵣₘₛ")
    plot!(w_rms[:],model.grid.zᵃᵃᶠ[1:Nz],label = "wᵣₘₛ")
    
end

gif(anim2, "rms_velocity.mp4", fps = 8)  # https://docs.juliaplots.org/latest/animations/

u_end=file["timeseries/u/$(iterations[end])"]
v_end=file["timeseries/v/$(iterations[end])"]
w_end=file["timeseries/w/$(iterations[end])"]

u_avg=mean(u_end,dims=(1,2));
v_avg=mean(v_end,dims=(1,2));
w_avg=mean(w_end,dims=(1,2));

u_prime=zeros(Nx,Ny,Nz)
v_prime=zeros(Nx,Ny,Nz)
w_prime=zeros(Nx,Ny,Nz)  #its dimension is actually Nx Ny Nz+1
for k in 1:Nz
    u_prime[:,:,k]=u_end[:,:,k].-u_avg[k]  #broadcast array.-scalar 
    v_prime[:,:,k]=v_end[:,:,k].-v_avg[k]
    w_prime[:,:,k]=w_end[:,:,k].-w_avg[k]
end

u_rms=(mean((u_end-u_prime).^2,dims=(1,2))).^0.5
v_rms=(mean((v_end-v_prime).^2,dims=(1,2))).^0.5
w_rms=(mean((w_end[:,:,1:Nz]-v_prime).^2,dims=(1,2))).^0.5;

plot(u_rms[:],model.grid.zᵃᵃᶠ[1:Nz],label = "uᵣₘₛ",
     ylabel = "Deepth (m)",
     xlabel = "Velocities (ms⁻¹)")
plot!(v_rms[:],model.grid.zᵃᵃᶠ[1:Nz],label = "vᵣₘₛ")
plot!(w_rms[:],model.grid.zᵃᵃᶠ[1:Nz],label = "wᵣₘₛ")

#close(file)

#=
using StatsBase
u_end64=convert(Array{Float64, 3}, u_end)
msd(u_end64,u_prime)
=#

keys(file)

model.velocities.w



grid.Δzᵃᵃᶜ;

file["grid"]["Δzᵃᵃᶜ"];






