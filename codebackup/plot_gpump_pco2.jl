
using OceanBioME, Oceananigans, Printf
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using NetCDF
using Statistics
using Plots
using Interpolations

Lx=20
Ly=20
Nz=50
########################################################## temperature
t_temperature_node=[0.,66,95,240,364]
temperature_idealize=[9.0,8.05,8.05,13.65,9.0]
temperature_itp = LinearInterpolation((t_temperature_node)days, temperature_idealize) 
t_function(x, y, z, t) = temperature_itp(mod(t, 364days)) 
########################################################## 
########################################################## salinity
s_function(x, y, z, t) = 35.
########################################################## MLD
t_mldplus_node=[0.,55,85,100,300,364]
mldplus_idealize=[280,420,420,40,40,280]
##########################################################  Mercator_physics
Mercator_physics = "subpolar_Si.nc" #A small sample of data downloaded is stored in subpolar.nc for ease of use. this is the forcing data 
time = ncread(Mercator_physics, "time")    # time in seconds
mld = ncread(Mercator_physics, "mld")      # mixed layer depth in Meters


########################################################## pco2
function pco2(result,output)
    for (iter, t) in enumerate(result.t)
        DIC = result.results[10, 1, 1, Nz, iter]   #["NO₃", "NH₄", "P", "Z", "D", "DD", "Dᶜ", "DDᶜ", "DOM", "DIC", "ALK", "OXY", "PAR"]
        ALK = result.results[11, 1, 1, Nz, iter]
        T = t_function(0.5*Lx, 0.5*Ly, 0.0, t) + 273.15   # here there is no 273.15 for flux , but has 273.15 for pco2 calculation 
        S = s_function(0.5*Lx, 0.5*Ly, 0.0, t)
        push!(output, Boundaries.pCO₂(DIC, ALK, T, S, Boundaries.defaults.airseaflux))
        # push!(output, Boundaries.airseaflux(0.5*Lx,0.5*Ly, t, DIC, ALK, T, S, Boundaries.defaults.airseaflux)) #function pCO₂(DIC, ALK, T, S, params)
    end
end
########################################################## CO2 flux

function flux_co2(result,output)                         # mmolC/m²s 
    for (iter, t) in enumerate(result.t)
        DIC = result.results[10, 1, 1, Nz, iter]   #tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DIC,:ALK)
        ALK = result.results[11, 1, 1, Nz, iter]
        T = t_function(0.5*Lx, 0.5*Ly, 0.0, t)   # here there is no 273.15
        S = s_function(0.5*Lx, 0.5*Ly, 0.0, t)
        gas = :CO₂
        push!(output, Boundaries.airseaflux(NaN, NaN, NaN, DIC, ALK, T, S, merge(Boundaries.defaults.airseaflux, (conc_air = (CO₂ = 400,), ),(;gas ))))  # merge(Boundaries.defaults.airseaflux, (conc_air = (CO₂ = 400,), ))
        # push!(output, Boundaries.airseaflux(0.5*Lx,0.5*Ly, t, DIC, ALK, T, S, Boundaries.defaults.airseaflux)) #function pCO₂(DIC, ALK, T, S, params)
    end
end

########################################################## NPP
function NPP(result,output)     
    for (iter, t) in enumerate(result.t)
        ##["NO₃", "NH₄", "P", "Z", "D", "DD", "Dᶜ", "DDᶜ", "DOM", "DIC", "ALK", "OXY", "PAR"]
        for j in 1:Nz
            NO₃ = result.results[1, 1, 1, j, iter] 
            NH₄ = result.results[2, 1, 1, j, iter] 
            P = result.results[3, 1, 1, j, iter] 
            Z = result.results[4, 1, 1, j, iter] 
            D = result.results[5, 1, 1, j, iter] 
            DD = result.results[6, 1, 1, j, iter] 
            Dᶜ = result.results[7, 1, 1, j, iter]   
            DDᶜ = result.results[8, 1, 1, j, iter]
            DOM = result.results[9, 1, 1, j, iter]
            PAR = result.results[13, 1, 1, j, iter]
            #NPP_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, params)
            output[j,iter] = LOBSTER.NPP_forcing(0.5*Lx, 0.5*Ly, result.z[j], t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, LOBSTER.defaults)
        end
    end
end



include("old_plotting_uitils.jl")

fs=8
xlabel="Time (days)"
title_pco2="pco2 (ppm)"

# ########################################################## Mercator BGC data 
# Mercator_BGC = "/store/DAMTP/sc2350/BGC.jl/subpolar/copernicus/subpolar_BGC.nc"  #Global ocean biogeochemistry forecast    2020-1-1: 2022-11-1 /store/DAMTP/sc2350/copernicus/

# datalength = 731 #  1036 #days
# ################################################### spco2 Pa
# spco2 = ncread(Mercator_BGC, "spco2");  # in Pa  
# spco2_mean = mean(spco2, dims=(1,2))[1,1,:] 
# spco2_std = std(spco2, dims=(1,2))[1,1,:] 
# Mercator_spco2 = Plots.plot(spco2_mean,ribbon=spco2_std,xlims=(0,730),label = "Mercator",xlabel=xlabel,ylabel="pCO2 (Pa)",legend=:bottomright)
# # savefig("pco2_forecastdata.pdf")

##########################################################  simulation 
path0="./withoutkelp/"
path1="./kelp01/"
path2="./kelp1/"
path3="./kelp10/"
path4="./kelp10_005/"

path5="./continuous_density10_largescale_couple_new/"
file5= path5*"withkelp_density10_largescale_couple_new.nc"   # 
DIC5=ncread(file5, "DIC")   

file0 = jldopen(path0*"withoutkelp.jld2", "r")
file1 = jldopen(path1*"kelp_01.jld2", "r")
file2 = jldopen(path2*"kelp_1.jld2", "r")
file3 = jldopen(path3*"kelp_10.jld2", "r")
file4 = jldopen(path4*"kelp_10.jld2", "r")

results0 = load_tracers(file0)         
datalength0 = length(results0.results[1,1,1,1,:])
results1 = load_tracers(file1)
results2 = load_tracers(file2)
results3 = load_tracers(file3)
results4 = load_tracers(file4)
pco2_0=[]
pco2_1=[]
pco2_2=[]
pco2_3=[]
pco2_4=[]

pco2(results0,pco2_0)
pco2(results1,pco2_1)
pco2(results2,pco2_2)
pco2(results3,pco2_3)
pco2(results4,pco2_4)

flux_0=[]
flux_1=[]
flux_2=[]
flux_3=[]
flux_4=[]
flux_co2(results0,flux_0)
flux_co2(results1,flux_1)
flux_co2(results2,flux_2)
flux_co2(results3,flux_3)
flux_co2(results4,flux_4)

Plots.plot(flux_0*Lx*Ly)  #mmolC/m²s  to mmolC/s 
Plots.plot!(flux_1*Lx*Ly)
Plots.plot!(flux_2*Lx*Ly)
Plots.plot!(flux_3*Lx*Ly)
Plots.plot!(flux_4*Lx*Ly)

Plots.plot(pco2_0)
Plots.plot!(pco2_1)
Plots.plot!(pco2_2)
Plots.plot!(pco2_3)
Plots.plot!(pco2_4)

Plots.plot(pco2_0, xlabel=xlabel,ylabel=title_pco2,label = "without kelp",legend=:bottomleft,xlims=(0,2190),ylims=(0,500),color="gray",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)
Plots.plot!( pco2_1, label = "kelp density 0.1",color="green")
Plots.plot!( pco2_2, label = "kelp density 1",color="blue")
Plots.plot!( pco2_3, label = "kelp density 10",color="red")
Plots.plot!( pco2_4, label = "kelp density 10, c=0.05",color="black")
# savefig("pco2.pdf")


############################################ Gravitational pump (mgC m^-2 d^-1)

ylabel_gpump="Gravitational pump (mgC m^-2 d^-1)"
title_gpump="at 100m depth"

V_d = -3.47e-5  #  Detritus sedimentation speed   ms⁻¹
V_dd = -200/day  #  Detritus sedimentation speed  -v_dd_min       50m/day=0.0005878  ms⁻¹ 200/day=2.3e-3m/s
molar_mass_C=12 # 12g/mol
z_index = 42  # julia> results0.z[9] -498.0,   julia> results0.z[42] -102.0

#["NO₃", "NH₄", "P", "Z", "D", "DD", "Dᶜ", "DDᶜ", "DOM", "DIC", "ALK", "OXY", "PAR"]

Dc_0 = results0.results[7, 1, 1, z_index, :]      
DDc_0 = results0.results[8, 1, 1, z_index, :]  
Dc_1 = results1.results[7, 1, 1, z_index, :] 
DDc_1 = results1.results[8, 1, 1, z_index, :] 
Dc_2 = results2.results[7, 1, 1, z_index, :] 
DDc_2 = results2.results[8, 1, 1, z_index, :] 
Dc_3 = results3.results[7, 1, 1, z_index, :] 
DDc_3 = results3.results[8, 1, 1, z_index, :] 
Dc_4 = results4.results[7, 1, 1, z_index, :] 
DDc_4 = results4.results[8, 1, 1, z_index, :] 

G_pump0 = (Dc_0*V_d+DDc_0*V_dd)*molar_mass_C*1day
G_pump1 = (Dc_1*V_d+DDc_1*V_dd)*molar_mass_C*1day
G_pump2 = (Dc_2*V_d+DDc_2*V_dd)*molar_mass_C*1day
G_pump3 = (Dc_3*V_d+DDc_3*V_dd)*molar_mass_C*1day
G_pump4 = (Dc_4*V_d+DDc_4*V_dd)*molar_mass_C*1day


Plots.plot(G_pump0, xlabel=xlabel,ylabel=ylabel_gpump,title=title_gpump,label = "without kelp",legend=:bottomleft,xlims=(0,2190),ylims=(-1000,0),color="gray",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)
Plots.plot!(G_pump1,label = "kelp density 0.1",color="green")
Plots.plot!(G_pump2, label = "kelp density 1",color="blue")
Plots.plot!(G_pump3, label = "kelp density 10",color="red")
Plots.plot!(G_pump4, label = "kelp density 10, c=0.05",color="black")


# savefig("G_pump@100m.pdf")


# p4=load_particles(path4*"kelp_10_1_particles.jld2")
# ww=particles(p4)
# ww[1]