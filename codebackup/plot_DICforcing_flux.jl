
using OceanBioME, Oceananigans, Printf
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using NetCDF
using Statistics
using Plots
using Interpolations

Lx=20
Ly=20
Nz=50

dx=20.0
dy=20.0
dz=12.0
vcell = dx*dy*dz
acell = dx*dy 
duration=6years
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
##########################################################  DIC_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, DIC, ALK, params)

function DICforcing(result,output)
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
            DIC = result.results[10, 1, 1, j, iter]
            ALK = result.results[11, 1, 1, j, iter]
            # OXY = result.results[12, 1, 1, j, iter]
            PAR = result.results[13, 1, 1, j, iter]
            # output[j,floor(Int,t/days)+1] = LOBSTER.DIC_forcing(0.5*Lx, 0.5*Ly, result.z[j], t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR, LOBSTER.defaults)
            output[j,iter] = LOBSTER.DIC_forcing(0.5*Lx, 0.5*Ly, result.z[j], t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, DIC, ALK, LOBSTER.defaults)
        end
    end
end



include("old_plotting_uitils.jl")

fs=8
xlabel="Time (days)"
ylabel_DICforcing="integrated DIC forcing (mmolC/s)"
ylabel_flux="integrated flux (mmolC/s)"
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

# path5="./continuous_density10_largescale_couple_new/"
# file5= path5*"withkelp_density10_largescale_couple_new.nc"   # 
# DIC5=ncread(file5, "DIC")   

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

fluxsum0 = flux_0*acell
fluxsum1 = flux_1*acell
fluxsum2 = flux_2*acell
fluxsum3 = flux_3*acell
fluxsum4 = flux_4*acell


Plots.plot(-fluxsum0,xlabel=xlabel,ylabel=ylabel_flux,label = "integrated flux without kelp",legend=:bottomleft,xlims=(0,2190),ylims=(-0.5,0.2),line=(:dot,:gray),titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #mmolC/m²s  to mmolC/s 
Plots.plot!(-fluxsum3,label = "integrated flux kelp density 10",line=(:dot,:red))
# savefig("integrated_flux.pdf")


############################################ DICforcing (mmolC/s)

DICforcing_0 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results0,DICforcing_0)
DICforcing_1 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results1,DICforcing_1)
DICforcing_2 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results2,DICforcing_2)
DICforcing_3 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results3,DICforcing_3)
DICforcing_4 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results4,DICforcing_4)
# ,right_margin=20Plots.mm

depth_index=1:Nz   #9:Nz 
DICforcingsum0 = sum((DICforcing_0.*vcell)[depth_index,:],dims=1)[:]  
DICforcingsum1 = sum((DICforcing_1.*vcell)[depth_index,:],dims=1)[:]  
DICforcingsum2 = sum((DICforcing_2.*vcell)[depth_index,:],dims=1)[:]  
DICforcingsum3 = sum((DICforcing_3.*vcell)[depth_index,:],dims=1)[1:1387]  
DICforcingsum4 = sum((DICforcing_4.*vcell)[depth_index,:],dims=1)[:]  

Plots.plot!(DICforcingsum0,xlabel=xlabel,ylabel=ylabel_DICforcing,label = "integrated DIC forcing without kelp",legend=:bottomleft,xlims=(0,2190),ylims=(-0.5,0.2),color="gray",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)
Plots.plot!(DICforcingsum3,label = "integrated DIC forcing kelp density 10",color="red")
# savefig("integrated_DICforcing_flux.pdf")






# p4=load_particles(path4*"kelp_10_1_particles.jld2")
# ww=particles(p4)
# ww[1]