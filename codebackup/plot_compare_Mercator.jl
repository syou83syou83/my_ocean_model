
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
 

########################################################## Mercator BGC data 
Mercator_BGC = "/store/DAMTP/sc2350/BGC.jl/subpolar/copernicus/subpolar_BGC.nc"  #Global ocean biogeochemistry forecast    2020-1-1: 2022-11-1 /store/DAMTP/sc2350/copernicus/

datalength = 731 #  1036 #days
#################################################### P    mmol m-3 
Rd_phy = 6.56 # C:N ratio for P molC molN-1
V_d = -3.47e-5  #  Detritus sedimentation speed   ms⁻¹
V_dd = -200/day  #  Detritus sedimentation speed  -v_dd_min       50m/day=0.0005878  ms⁻¹ 200/day=2.3e-3m/s
molar_mass_C=12 # 12g/mol

phyc = ncread(Mercator_BGC, "phyc");  #phyc
phy_mean = mean(phyc, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
depth_phyc = ncread(Mercator_BGC, "depth");
Mercator_P = Plots.heatmap(1:datalength, -depth_phyc[end:-1:1], phy_mean[end:-1:1,1:datalength],xlims=(0,730),ylims=(-600,0),xlabel=xlabel,ylabel=ylabel,title=title_P)  #1:730
Mercator_P=Plots.plot!(0:729,-[mld;mld],xlims=(0,730),ylims=(-600,0),label = "MLD_Mercator",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,legend=:bottomleft)
Mercator_P=Plots.plot!([t_mldplus_node;t_mldplus_node.+365],-[mldplus_idealize;mldplus_idealize],xlims=(0,730),ylims=(-600,0),label = "MLD_Idealized",legend=:bottomright)
#################################################### NPP   mg m-3 day-1 (milligrams of Carbon per cubic meter per day )
nppv = ncread(Mercator_BGC, "nppv");  
nppv_mean = mean(nppv, dims=(1,2))[1,1,:,:] #
Mercator_NPP = Plots.heatmap(1:datalength, -depth_phyc[end:-1:1], nppv_mean[end:-1:1,1:datalength],ylims=(-600,0),titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_NPP)  #1:730

#################################################### N  mmol m-3 
no3 = ncread(Mercator_BGC, "no3");  #
no3_mean = mean(no3, dims=(1,2))[1,1,:,:] #
Mercator_N = Plots.heatmap(1:datalength, -depth_phyc[end:-1:1], no3_mean[end:-1:1,1:datalength],ylims=(-600,0),titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_NO₃)  #1:730
################################################### spco2 Pa
spco2 = ncread(Mercator_BGC, "spco2");  # in Pa  
spco2_mean = mean(spco2, dims=(1,2))[1,1,:] 
spco2_std = std(spco2, dims=(1,2))[1,1,:] 
Mercator_spco2 = Plots.plot(spco2_mean,ribbon=spco2_std,xlims=(0,730),label = "Mercator",xlabel=xlabel,ylabel="pCO2 (Pa)",legend=:bottomright)
# savefig("pco2_forecastdata.pdf")


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

Mercator_spco2 = Plots.plot!(1:datalength, pco2_0[1:datalength]/10,titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel="pCO2 (Pa)",label = "Simulation")
# savefig("subpolar_spco2.pdf")

flux_0=[]
flux_co2(results0,flux_0)
Plots.plot(flux_0*Lx*Ly)  #mmolC/m²s  to mmolC/s 

Plots.plot!(pco2_1/10)
Plots.plot!(pco2_2/10)
Plots.plot!(pco2_3/10)
Plots.plot!(pco2_4/10)

###########

####################################################### P simulation 
#["NO₃", "NH₄", "P", "Z", "D", "DD", "Dᶜ", "DDᶜ", "DOM", "DIC", "ALK", "OXY", "PAR"]
P_0 = results0.results[3, 1, 1, :, :]   
Sim_P = Plots.heatmap(1:datalength, results0.z, P_0[:,1:datalength],xlims=(0,730),ylims=(-600,0),titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_P)  #1:730
Sim_P = Plots.plot!([t_mldplus_node;t_mldplus_node.+365],-[mldplus_idealize;mldplus_idealize],label = "MLD_Idealized",legend=:bottomright)

Plots.plot(Mercator_P,Sim_P)
# savefig("subpolar_P.pdf")
########### nppv simulation 
# NPP_0 = zeros(Nz,floor(Int,datalength0)+1)  
# NPP(results0,NPP_0)
# NPP_0=NPP_0*Rd_phy*molar_mass_C*1days   #  from mmolN/m^3/s         to    mg m-3 day-1 (milligrams of Carbon per cubic meter per day )
# Sim_NPP0 = Plots.heatmap(1:datalength, results0.z,NPP_0[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_NPP)

# datalength = 731

NPP_0 = zeros(Nz,floor(Int,datalength0)+1)  
NPP_3 = zeros(Nz,floor(Int,datalength3)+1)  
NPP(results0,NPP_0)
NPP(results3,NPP_3)
NPP_0=NPP_0*Rd_phy*molar_mass_C*1days   #  from mmolN/m^3/s         to    mg m-3 day-1 (milligrams of Carbon per cubic meter per day )
NPP_3=NPP_3*Rd_phy*molar_mass_C*1days 
Sim_NPP0 = Plots.heatmap(1:datalength0, results0.z,NPP_0[:,1:datalength0],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title="NPP without kelp (mgC m^-3 day^-1)")
Sim_NPP3 = Plots.heatmap(1:datalength3, results3.z,NPP_3[:,1:datalength3],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title="NPP kelp density=10 (mgC m^-3 day^-1)")
Plots.plot(Sim_NPP0,Sim_NPP3)




# Plots.plot(Mercator_NPP,Sim_NPP0)
# savefig("subpolar_NPP.pdf")

########### N simulation 
NO₃_0 = results0.results[1, 1, 1, :, :]   
Sim_NO₃ = Plots.heatmap(1:datalength, results0.z, NO₃_0[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_NO₃)  #1:730
Plots.plot(Mercator_N,Sim_NO₃)
# savefig("subpolar_NO3.pdf")

NH₄_0 = results0.results[2, 1, 1, :, :]   
Sim_NH₄ = Plots.heatmap(1:datalength, results0.z, NH₄_0[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_NH₄)  #1:730

Z_0 = results0.results[4, 1, 1, :, :] 
D_0 = results0.results[5, 1, 1, :, :] 
DD_0 = results0.results[6, 1, 1, :, :] 
# Dc_0 = results0.results[7, 1, 1, :, :] 
# DDc_0 = results0.results[8, 1, 1, :, :] 
DOM_0 = results0.results[9, 1, 1, :, :] 
DIC_0 = results0.results[10, 1, 1, :, :] 
ALK_0 = results0.results[11, 1, 1, :, :] 


Sim_NH₄ = Plots.heatmap(1:datalength, results0.z, NH₄_0[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_NH₄)  #1:730
Sim_Z = Plots.heatmap(1:datalength, results0.z, Z_0[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_Z)  #1:730
Sim_D = Plots.heatmap(1:datalength, results0.z, D_0[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_D)  #1:730
Sim_DD = Plots.heatmap(1:datalength, results0.z, DD_0[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_DD)  #1:730
Sim_DOM = Plots.heatmap(1:datalength, results0.z, DOM_0[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_DOM)  #1:730
Sim_DIC = Plots.heatmap(1:datalength, results0.z, DIC_0[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_DIC)  #1:730
Sim_ALK = Plots.heatmap(1:datalength, results0.z, ALK_0[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_ALK)  #1:730

Plots.plot(Sim_NO₃, Sim_NH₄, Sim_P, Sim_Z, Sim_D, Sim_DD,Sim_DOM, Sim_DIC, Sim_ALK)
# savefig("tracers_without_kelp.pdf")


NO₃_3 = results3.results[1, 1, 1, :, :]   
NH₄_3 = results3.results[2, 1, 1, :, :]   
P_3 = results3.results[3, 1, 1, :, :]   
Z_3 = results3.results[4, 1, 1, :, :] 
D_3 = results3.results[5, 1, 1, :, :] 
DD_3 = results3.results[6, 1, 1, :, :] 
DOM_3 = results3.results[9, 1, 1, :, :] 
DIC_3 = results3.results[10, 1, 1, :, :] 
ALK_3 = results3.results[11, 1, 1, :, :] 

Sim_NO₃_3 = Plots.heatmap(1:datalength, results0.z, NO₃_3[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_NO₃)  #1:730
Sim_NH₄_3 = Plots.heatmap(1:datalength, results0.z, NH₄_3[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_NH₄)  #1:730
Sim_P_3 = Plots.heatmap(1:datalength, results0.z, P_3[:,1:datalength],xlims=(0,730),ylims=(-600,0),titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_P)  #1:730
Sim_Z_3 = Plots.heatmap(1:datalength, results0.z, Z_3[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_Z)  #1:730
Sim_D_3 = Plots.heatmap(1:datalength, results0.z, D_3[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_D)  #1:730
Sim_DD_3 = Plots.heatmap(1:datalength, results0.z, DD_3[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_DD)  #1:730
Sim_DOM_3 = Plots.heatmap(1:datalength, results0.z, DOM_3[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_DOM)  #1:730
Sim_DIC_3 = Plots.heatmap(1:datalength, results0.z, DIC_3[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_DIC)  #1:730
Sim_ALK_3 = Plots.heatmap(1:datalength, results0.z, ALK_3[:,1:datalength],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,title=title_ALK)  #1:730
Plots.plot(Sim_NO₃_3, Sim_NH₄_3, Sim_P_3, Sim_Z_3, Sim_D_3, Sim_DD_3, Sim_DOM_3, Sim_DIC_3, Sim_ALK_3)
# savefig("tracers_with_kelp.pdf")
