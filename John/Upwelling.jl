# A test of a simple upwelling system

using Printf
using Oceananigans
using Oceananigans.Units: hours, hour, day, days, minute, minutes, second, seconds

# ## The grid

Lx=200e3;
Ly=200e3;
Lz=2e3;
Nx=128;
Ny=128;
Nz=32;
duration=120days;

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
                                         z = (0, Lz),
#                                         z = z_faces,
                                          topology = (Bounded, Periodic, Bounded))

# ## Rotation
#
# The classical Eady problem is posed on an ``f``-plane. We use a Coriolis parameter
# typical to mid-latitudes on Earth,

coriolis = FPlane(f=-1e-4) # [s⁻¹]

# ## Buoyancy that depends on temperature and salinity
#
# We use the `SeawaterBuoyancy` model with a linear equation of state,

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=2e-4, β=0))


# ## The background flow
#
# We build a `NamedTuple` of parameters that describe the background flow,

basic_state_parameters = ( f = coriolis.f,      # s⁻¹, Coriolis parameter
                           N = 1e-4,            # s⁻¹, buoyancy frequency
                           Lz = grid.Lz)         # m, ocean depth

# and then construct the background fields ``U`` and ``B``

# ## Boundary conditions
#

# u₁₀ = 10    # m s⁻¹, average wind velocity 10 meters above the ocean
# cᴰ = 2.5e-3 # dimensionless drag coefficient
# ρₐ = 1.225  # kg m⁻³, average density of air at sea-level
# cᴾ = 3991 # J K⁻¹ kg⁻¹, typical heat capacity for seawater
#Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²

ρₒ = 1026 # kg m⁻³, reference density
τʷ=0.1
Qᵘ=-τʷ/ρₒ

cᴰ = 6e-4 # linear drag coefficient
@inline drag_u(x, y, t, u, v, cᴰ) = - cᴰ * u 
@inline drag_v(x, y, t, u, v, cᴰ) = - cᴰ * v 
drag_bc_u = FluxBoundaryCondition(drag_u, field_dependencies=(:u, :v), parameters=cᴰ)
drag_bc_v = FluxBoundaryCondition(drag_v, field_dependencies=(:u, :v), parameters=cᴰ)
u_bcs = FieldBoundaryConditions(bottom = drag_bc_u)
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ), bottom = drag_bc_v)

# Configure Lagrangian particles
n_particles = 100;

x₀ = 125*rand(n_particles);
y₀ = 125*rand(n_particles);
z₀ = 0*ones(n_particles);

lagrangian_particles = LagrangianParticles(x=x₀, y=y₀, z=z₀, restitution=0)

# Tracer forcing
#C_forcing_function(x, y, z, t, C, params) = -params.m * C + params.μ₀*exp((-(x-Lx*9/10)^2-(y-Ly/2)^2)/(1e3)^2)*(1-C)
C_forcing_function(x, y, z, t, C, params) = -params.m * C + 0.5*params.μ₀*(1+tanh((x-(Lx-2e4))/1e4))*(1-C)

#C_forcing_function(x, y, z, t, C, params) = params.μ₀*(0.5*(1+tanh((x-(Lx-2e4))/1e4))-C)
nothing # hide
# with parameters
C_forcing_parameters = (μ₀ = 10/day,   # surface growth rate
                                 m = 0.1/day) # mortality rate due to virus and zooplankton grazing
# We tell `Forcing` that our forcing depends
# on the concentration `c` and the chosen parameters,
C_forcing = Forcing(C_forcing_function, field_dependencies = :C,
                            parameters = C_forcing_parameters)

# Immersed boundary forcing
# Define the bottom depth
xₛ=100e3; # location of mid-shelf break
Lₛ=50e3; # width of continental slope
@inline bottom(x,y)=-Lz*(1-tanh((x-xₛ)/Lₛ))/2;
# Define the mask function
@inline mask(x,y,z,params)=(1-tanh((z-bottom(x,y))/params.mask_width))/2;

@inline IB_u_function(x, y, z, t, u, params) = @inbounds -mask(x, y, z, params) * u / params.tau 
@inline IB_v_function(x, y, z, t, v, params) = @inbounds -mask(x, y, z, params) * v / params.tau 
@inline IB_w_function(x, y, z, t, w, params) = @inbounds -mask(x, y, z, params) * w / params.tau 

IB_parameters = (mask_width=20, # vertical thickness of transition at the edge of the mask
                tau=30minutes) # forcing timescale
IB_u_forcing = Forcing(IB_u_function, field_dependencies = :u, parameters=IB_parameters)
IB_v_forcing = Forcing(IB_v_function, field_dependencies = :v, parameters=IB_parameters)
IB_w_forcing = Forcing(IB_w_function, field_dependencies = :w, parameters=IB_parameters)

# ## Turbulence closures
#
# We use a horizontal hyperdiffusivity and a Laplacian vertical diffusivity
# to dissipate energy in the Eady problem.
# To use both of these closures at the same time, we set the keyword argument
# `closure` to a tuple of two closures.

κ₂z = 1e-5 # [m² s⁻¹] Laplacian vertical viscosity and diffusivity
κ₄h = 1e-4 / day * (Lx/Nx)^4 # [m⁴ s⁻¹] horizontal hyperviscosity and hyperdiffusivity

Laplacian_vertical_diffusivity = AnisotropicDiffusivity(νh=0, κh=0, νz=κ₂z, κz=κ₂z)
biharmonic_horizontal_diffusivity = AnisotropicBiharmonicDiffusivity(νh=κ₄h, κh=κ₄h)

# ## Model instantiation
#
# We instantiate the model with the fifth-order WENO advection scheme, a 3rd order
# Runge-Kutta time-stepping scheme, and a `BuoyancyTracer`.

model = NonhydrostaticModel(
           #architecture = CPU(),
                   grid = grid,
              advection = UpwindBiasedFifthOrder(),
            timestepper = :RungeKutta3,
               coriolis = coriolis,
                tracers = (:T, :S, :C),
               buoyancy = buoyancy,
               particles = lagrangian_particles,
                closure = (Laplacian_vertical_diffusivity, biharmonic_horizontal_diffusivity),
                forcing = (C=C_forcing,u=IB_u_forcing, v=IB_v_forcing, w=IB_w_forcing),
    boundary_conditions = (u=u_bcs, v=v_bcs)
)

# ## Initial conditions
#
# We seed our initial conditions with random noise stimulate the growth of
# baroclinic instability.

## A noise function, damped at the top and bottom
Ξ(z) = randn();
## Velocity initial condition: random noise scaled by the friction velocity.
uᵢ(x, y, z) = 1e-1 * Ξ(z);

Tₐ=4;
Tᵦ=15.5;
h₁=80;
h₂=-80;
dTdx=-5/600e3;
## Temperature initial condition
Tᵢ(x, y, z) = Tₐ*(1+tanh((z-h₂)/h₁))+Tᵦ*exp(z/1000)+dTdx*x

vᵢ(x, y, z) = 2e-4*9.81*dTdx*(z+Lz)/-1e-4;

dSdz=1/1000
Sᵢ(x, y, z) = 35.5 + dSdz * z

Cᵢ(x, y, z) = 0.0

set!(model, u=uᵢ, v=vᵢ, T=Tᵢ, S=Sᵢ, C=Cᵢ)

# We subtract off any residual mean velocity to avoid exciting domain-scale
# inertial oscillations. We use a `sum` over the entire `parent` arrays or data
# to ensure this operation is efficient on the GPU (set `architecture = GPU()`
# in `NonhydrostaticModel` constructor to run this problem on the GPU if one
# is available).

ū = sum(model.velocities.u.data.parent) / (grid.Nx * grid.Ny * grid.Nz)
v̄ = sum(model.velocities.v.data.parent) / (grid.Nx * grid.Ny * grid.Nz)

model.velocities.u.data.parent .-= ū
model.velocities.v.data.parent .-= v̄
nothing # hide

# ## Simulation set-up
#
simulation = Simulation(model, Δt=10minutes, stop_time=duration)

# The `TimeStepWizard` helps ensure stable time-stepping
# with a Courant-Freidrichs-Lewy (CFL) number of 1.0.

wizard = TimeStepWizard(cfl=1, max_change=1.1, max_Δt=10minutes)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# ### A progress messenger
#
# We add a callback that prints out a helpful progress message while the simulation runs.

start_time = time_ns()

progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(1e-9 * (time_ns() - start_time)),
                        prettytime(sim.Δt),
                        AdvectiveCFL(sim.Δt)(sim.model))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# ### Output
#
# To visualize the baroclinic turbulence ensuing in the Eady problem,
# we use `ComputedField`s to diagnose and output vertical vorticity and divergence.
# Note that `ComputedField`s take "AbstractOperations" on `Field`s as input:

#u, v, w = model.velocities # unpack velocity `Field`s
## Vertical vorticity [s⁻¹]
#vorticity = ComputedField(∂x(v) - ∂y(u))
# Create a NamedTuple for vorticity
#vorticity_tuple = (; ω = vorticity)

simulation.output_writers[:xz_slices] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                          prefix = "slice_xz",
                     field_slicer = FieldSlicer(j=Int(grid.Ny/2)),
                         schedule = TimeInterval(1hours),
                            force = true)

 simulation.output_writers[:xy_slices] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                          prefix = "slice_xy",
                     field_slicer = FieldSlicer(k=Int(grid.Nz)),
                         schedule = TimeInterval(1hours),
                            force = true)

#  simulation.output_writers[:particles] = JLD2OutputWriter(model, (particles=model.particles,), 
#                             prefix = "particles",
#                           schedule = TimeInterval(1minute),
#                              force = true)
nothing # hide

# Define a checkpoint file to restart the model later
simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=TimeInterval(10days), prefix="model_checkpoint")

# All that's left is to press the big red button:

run!(simulation)

