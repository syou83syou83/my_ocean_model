
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
        push!(output, Boundaries.airseaflux(NaN, NaN, NaN, DIC, ALK, T, S, merge(Boundaries.defaults.airseaflux, (;gas ))))
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

fs=5
xlabel="Time (days)"
ylabel="z (m)"
title_logP="log(P)"
title_P="P (mmolN/m^3)"
title_NPP="NPP (mgC m^-3 day^-1)"
title_NO₃="NO₃ (mmolN/m^3)"
title_NH₄="NH₄ (mmolN/m^3)"
title_Z="Z (mmolN/m^3)"
title_D="D (mmolN/m^3)"
title_DD="DD (mmolN/m^3)"
title_DIC="DIC (mmolC/m^3)"
title_ALK="ALK (mmolN/m^3)"
title_DOM="DOM (mmolN/m^3)"
ylabel_DICforcing="integrated DIC forcing (mmolC/s)"
ylabel_flux="integrated flux (mmolC/s)"
##########################################################  simulation 
path0="./withoutkelp/"
path1="./kelp01/"
path2="./kelp1/"
path3="./kelp10/"
path4="./kelp10_005/"

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
datalength3 = length(results3.results[1,1,1,1,:])

results4 = load_tracers(file4)



###########

########### nppv simulation 
Rd_phy = 6.56 # C:N ratio for P molC molN-1
molar_mass_C=12 # 12g/mol

NPP_0 = zeros(Nz,floor(Int,datalength0)+1)  
NPP_3 = zeros(Nz,floor(Int,datalength3)+1)  
NPP(results0,NPP_0)
NPP(results3,NPP_3)
NPP_0=NPP_0*Rd_phy*molar_mass_C*1days   #  from mmolN/m^3/s         to    mg m-3 day-1 (milligrams of Carbon per cubic meter per day )
NPP_3=NPP_3*Rd_phy*molar_mass_C*1days 
Sim_NPP0 = Plots.heatmap(1:datalength0, results0.z,NPP_0[:,1:datalength0],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title="NPP without kelp (mgC m^-3 day^-1)")
Sim_NPP3 = Plots.heatmap(1:datalength3, results3.z,NPP_3[:,1:datalength3],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title="NPP kelp density=10 (mgC m^-3 day^-1)")
Plots.plot(Sim_NPP0,Sim_NPP3)

# savefig("NPP_kelp.pdf")

Z_0 = results0.results[4, 1, 1, :, :] 
D_0 = results0.results[5, 1, 1, :, :] 
DD_0 = results0.results[6, 1, 1, :, :] 
Dc_0 = results0.results[7, 1, 1, :, :] 
DDc_0 = results0.results[8, 1, 1, :, :] 
DOM_0 = results0.results[9, 1, 1, :, :] 
DIC_0 = results0.results[10, 1, 1, :, :] 
ALK_0 = results0.results[11, 1, 1, :, :] 


NO₃_3 = results3.results[1, 1, 1, :, :]   
NH₄_3 = results3.results[2, 1, 1, :, :]   
P_3 = results3.results[3, 1, 1, :, :]   
Z_3 = results3.results[4, 1, 1, :, :] 
D_3 = results3.results[5, 1, 1, :, :] 
DD_3 = results3.results[6, 1, 1, :, :] 
DOM_3 = results3.results[9, 1, 1, :, :] 
DIC_3 = results3.results[10, 1, 1, :, :] 
ALK_3 = results3.results[11, 1, 1, :, :] 


depth_index=1:Nz 

DIC0sum=sum((DIC_0.*vcell)[depth_index,:],dims=1)[:]
DIC3sum=sum((DIC_3.*vcell)[depth_index,:],dims=1)[:]
Plots.plot(DIC0sum,ylims=(5.2e8,5.4e8),ylabel="DIC sum (mmol)", line=(:solid,:blue),label="integrated DIC without kelp",legendfontsize=5,legend=:bottomright,right_margin=30Plots.mm)
Plots.plot!(DIC3sum,ylims=(5.2e8,5.4e8),ylabel="DIC sum (mmol)", line=(:solid,:black),label="integrated DIC without kelp",legendfontsize=5,legend=:bottomright,right_margin=30Plots.mm)

flux_0=[]
flux_3=[]

flux_co2(results0,flux_0)
flux_co2(results3,flux_3)

fluxsum0 = flux_0*acell
fluxsum3 = flux_3*acell


Plots.plot(-fluxsum0*1days,xlabel=xlabel,ylabel=ylabel_flux,label = "integrated flux without kelp",legend=:bottomleft,xlims=(0,2190),ylims=(-0.5,0.2),line=(:dot,:gray),titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #mmolC/m²s  to mmolC/s 
Plots.plot!(-fluxsum3*1days,label = "integrated flux kelp density 10",line=(:dot,:red))

((-fluxsum3*1days)[end]-(fluxsum3*1days)[1])/(DIC3sum[1]-DIC3sum[end])
