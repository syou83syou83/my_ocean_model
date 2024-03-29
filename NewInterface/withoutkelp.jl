# # Simple active particle example
# In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and active particles modelling the growth of sugar kelp. This demonstraits:
# - How to setup OceanBioME's biogeochemical models
# - How to setup light attenuation
# - How to add biologically active particles which interact with the biodeochemical model
# - How to include optional tracer sets (carbonate chemistry and oxygen)
# - How to visulise results

# This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, Printf, GLMakie"
# ```

# ## Model setup
# We load the packages and choose the default LOBSTER parameter set
using OceanBioME, Oceananigans,Printf
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Interpolations
# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic)

##PAR⁰(x, y, t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((mod(t, 365days)-200days)/50days)^2))) .+ 2
##H(t, t₀, t₁) = ifelse(t₀<t<t₁, 1.0, 0.0)
##fmld1(t) = H.(t, 50days, 365days).*(1 ./(1 .+exp.(-(t-100days)/(5days)))).*(1 ./(1 .+exp.((t .-330days)./(25days))))
##MLD(t) = (-10 .-340 .*(1 .-fmld1(364.99999days).*exp.(-t/25days).-fmld1.(mod.(t, 365days))))
##κₜ(x, y, z, t) = 1e-2*max(1-(z+MLD(t)/2)^2/(MLD(t)/2)^2,0)+1e-4; 
##t_function(x, y, z, t) = 2.4*cos(t*2π/year + 50day) + 10
##s_function(x, y, z, t) = 35.0

########################################################## temperature
t_temperature_node=[0.,66,95,240,364]
temperature_idealize=[9.0,8.05,8.05,13.65,9.0]
temperature_itp = LinearInterpolation((t_temperature_node)days, temperature_idealize) 
t_function(x, y, z, t) = temperature_itp(mod(t, 364days)) 
########################################################## 
########################################################## salinity
s_function(x, y, z, t) = 35.0
########################################################## 
########################################################## PAR
t_par_node=[0.,30,120,200,330,364]
par_idealize=[3,3,90,90,3,3]
PAR_itp = LinearInterpolation((t_par_node)days, par_idealize)
PAR⁰(x, y, t) = PAR_itp(mod(t, 364days)) 
########################################################## 
########################################################## diffusivity 
t_mldplus_node=[0.,55,85,100,300,364]
mldplus_idealize=[280,420,420,40,40,280]
mld_itp = LinearInterpolation((t_mldplus_node)days, mldplus_idealize)  #in seconds 
κₜ(x, y, z, t) = 8e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-4;
##########################################################

# ## Grid and PAR field
# Define the grid and an extra Oceananigans field for the PAR to be stored in
##Lx, Ly = 20, 20
##grid = RectilinearGrid(size=(1, 1, 50), extent=(Lx, Ly, 200)) 
##PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  
duration= 6years
Lx, Ly, Lz = 20, 20, 600
Nx, Ny, Nz = 1, 1, 50
grid = RectilinearGrid(size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz)) 
PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  
# ## Kelp Particle setup
#@info "Setting up kelp particles"
#n = 100 # number of kelp fronds
#z₀ = [-100:-1;]*1.0 # depth of kelp fronds   z₀ = [-100:-1;]*1.0   [-21:5:-1;]*1.0
# kelp_particles = SLatissima.setup(;n, 
#                                   x₀ = Lx/2, y₀ = Ly/2, z₀, 
#                                   A₀ = 0.1, N₀ = 0.01, C₀ = 0.1, 
#                                   latitude = 57.5,
#                                   scalefactor = 10.0, 
#                                   T = t_function, S = s_function, urel = 0.2, 
#                                   optional_tracers = (:NH₄, :DIC, :bPOM, :bPOC, :O₂, :DON, :DOC))

# Specify the boundary conditions for DIC and O₂ based on the air-sea CO₂ and O₂ flux
CO₂_flux = GasExchange(; gas = :CO₂, temperature = t_function, salinity = s_function, air_concentration = 400)
O₂_flux = GasExchange(; gas = :O₂, temperature = t_function, salinity = s_function)
model = NonhydrostaticModel(; grid,
                              advection = UpwindBiasedThirdOrder(),  #WENO(;grid)
                              timestepper = :RungeKutta3,
                              closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                              biogeochemistry = LOBSTER(; grid,
                                                          surface_phytosynthetically_active_radiation = PAR⁰,
                                                          carbonates = true,
                                                          oxygen = true,
                                                          variable_redfield = true,
                                                          open_bottom = true),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),
                                                     O₂ = FieldBoundaryConditions(top = O₂_flux), ),
                              auxiliary_fields = (; PAR),
                            #   particles = kelp_particles
                              )

Pᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.002           #in mmolN m^-3
Zᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.008           #in mmolN m^-3
sPONᵢ(x, y, z)=0                                                #in mmolN m^-3    Dᵢ
bPONᵢ(x, y, z)=0                                               #in mmolN m^-3     DDᵢ
sPOCᵢ(x, y, z)=0                                                #in mmolC m^-3    Dᶜᵢ
bPOCᵢ(x, y, z)=0                                               #in mmolC m^-3     DDᶜᵢ
NO₃ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*6+11.4                #in mmolN m^-3
NH₄ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*0.05+0.05             #in mmolN m^-3
DONᵢ(x, y, z)= 0                                             #in mmolN m^-3
DOCᵢ(x, y, z)= 0                                             #in mmolC m^-3
DICᵢ(x, y, z)= 2200                                          #in mmolC m^-3
Alkᵢ(x, y, z)= 2400                                          #in mmolN m^-3
O₂ᵢ(x, y, z) = 240                                          #in mmolO m^-3

#set!(model, P=0.03, Z=0.03, NO₃=11.0, NH₄=0.05, DIC=2200.0, Alk=2400.0, O₂=240.0)
set!(model, P=Pᵢ, Z=Zᵢ, sPON=sPONᵢ, bPON=bPONᵢ, sPOC=sPOCᵢ, bPOC=bPOCᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DON=DONᵢ, DOC=DOCᵢ, DIC=DICᵢ, Alk=Alkᵢ, O₂=O₂ᵢ)  #没有速度了
#julia> keys(pp["timeseries"]) 15-element Vector{String}:  "NO₃"  "NH₄" "P" "Z" "sPON" "bPON" "DON" "DIC" "Alk" "O₂" "sPOC" "bPOC" "DOC" "PAR" "t"#

# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Couples the particles to the biodeochemical model
# - Update the PAR field from the surface PAR and phytoplankton concentration
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)

simulation = Simulation(model, Δt=2.5minutes, stop_time=duration) 

#simulation.callbacks[:couple_particles] = Callback(Particles.infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "withoutkelp"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true)
#simulation.output_writers[:particles] = JLD2OutputWriter(model, (particles=model.particles, ), filename = "$(filename)_particles.jld2", schedule = TimeInterval(1day), overwrite_existing = true)

simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON), warn=false), callsite = UpdateStateCallsite())
simulation.callbacks[:neg2] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON), warn=false), callsite = TendencyCallsite())
simulation.callbacks[:timestep] = Callback(update_timestep!, IterationInterval(1), (c_forcing=0.1, c_adv=0.5, c_diff=0.5, w = 200/day, relaxation=0.95), TimeStepCallsite())
simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=SpecifiedTimes([i*years for i=1:12]), prefix=filename) #prefix="kelp_01"
# ## Run!
# Finally we run the simulation
run!(simulation)

# Now we can visulise the results with some post processing to diagnose the air-sea CO₂ flux - hopefully this looks different to the example without kelp!

P = FieldTimeSeries("$filename.jld2", "P")
NO₃ = FieldTimeSeries("$filename.jld2", "NO₃")
Z = FieldTimeSeries("$filename.jld2", "Z")
sPON = FieldTimeSeries("$filename.jld2", "sPON") 
bPON = FieldTimeSeries("$filename.jld2", "bPON")
DIC = FieldTimeSeries("$filename.jld2", "DIC")
sPOC = FieldTimeSeries("$filename.jld2", "sPOC")
bPOC = FieldTimeSeries("$filename.jld2", "bPOC")
Alk = FieldTimeSeries("$filename.jld2", "Alk")

x, y, z = nodes(P)
times = P.times

air_sea_CO₂_flux = zeros(size(P)[4])
carbon_export = zeros(size(P)[4])
for (i, t) in enumerate(times)
    air_sea_CO₂_flux[i] = CO₂_flux.condition.parameters(0.0, 0.0, t, DIC[1, 1, 50, i], Alk[1, 1, 50, i], t_function(1, 1, 0, t), s_function(1, 1, 0, t))*Oceananigans.Operators.Ax(1, 1, 50, grid, Center(), Center(), Center())
    carbon_export[i] = (sPOC[1, 1, end-20, i]*model.biogeochemistry.sinking_velocities.sPOM.w[1] .+ bPOC[1, 1, end-20, i]*model.biogeochemistry.sinking_velocities.bPOM.w[1])
end

using GLMakie
f=Figure(backgroundcolor=RGBf(1, 1, 1), fontsize=30)

axP = Axis(f[1, 1:2], ylabel="z (m)", xlabel="Time (days)", title="Phytoplankton concentration (mmol N/m³)")
hmP = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(P[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbP = Colorbar(f[1, 3], hmP)

axNO₃ = Axis(f[1, 4:5], ylabel="z (m)", xlabel="Time (days)", title="Nitrate concentration (mmol N/m³)")
hmNO₃ = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(NO₃[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbNO₃ = Colorbar(f[1, 6], hmNO₃)

axZ = Axis(f[2, 1:2], ylabel="z (m)", xlabel="Time (days)", title="Zooplankton concentration (mmol N/m³)")
hmZ = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(Z[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbZ = Colorbar(f[2, 3], hmZ)

axD = Axis(f[2, 4:5], ylabel="z (m)", xlabel="Time (days)", title="Detritus concentration (mmol C/m³)")
hmD = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(sPOC[1, 1, end-23:end, 1:101])' .+ float.(bPOC[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbD = Colorbar(f[2, 6], hmD)

axfDIC = Axis(f[3, 1:4], xlabel="Time (days)", title="Air-sea CO₂ flux and Sinking", ylabel="Flux (kgCO₂/m²/year)")
hmfDIC = GLMakie.lines!(times./days, air_sea_CO₂_flux.*(12+16*2).*year/(1000*1000), label="Air-sea flux")
hmfExp = GLMakie.lines!(times./days, carbon_export.*(12+16*2).*year/(1000*1000), label="Sinking export")

f[3, 5] = Legend(f, axfDIC, "", framevisible = false)

save("$(filename).png", f)

