

using JLD2
# using Printf
# using Oceananigans
using Plots
# using Statistics
using Interpolations
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using NetCDF


path0="./withoutkelpbase_realforcing/"
file0= path0*"withoutkelp_long_realforcing.nc"   #    

path1="./withoutkelpbase/"
file1= path1*"withoutkelp_long.nc"   #    
path2="./continuous_density01_largescale_couple/"
file2= path2*"withkelp_density01_largescale_couple.nc"   # density=0.1
path3="./continuous_density1_largescale_couple/"
file3= path3*"withkelp_density1_largescale_couple.nc"   # 
path4="./continuous_density10_largescale_couple/"
file4= path4*"withkelp_density10_largescale_couple.nc"   # 
path5="./continuous_density100_largescale_couple/"
file5= path5*"withkelp_density100_largescale_couple.nc"   # 
path6="./pickup_density1000_largescale_couple/"
file6= path6*"withkelp_density1000_largescale.nc"   # 
path7="./pickup_density10000_largescale_couple/"
file7= path7*"withkelp_density10000_largescale.nc"   # 
path8="./pickup_density10_largescale_couple/"
file8= path8*"withkelp_density10_largescale.nc"   # 

z0=ncread(file0, "zC")
t0=ncread(file0, "time")/1days
DD0=ncread(file0, "DD")          
D0=ncread(file0, "D")
Phyc0=ncread(file0, "P")
NO₃0=ncread(file0, "NO₃")

z1=ncread(file1, "zC")
t1=ncread(file1, "time")/1days
DD1=ncread(file1, "DD")          #DD1[1,1,:,:]
D1=ncread(file1, "D")
Phyc1=ncread(file1, "P")
NO₃1=ncread(file1, "NO₃")

t2=ncread(file2, "time")/1days
DD2=ncread(file2, "DD")
D2=ncread(file2, "D")

t3=ncread(file3, "time")/1days
DD3=ncread(file3, "DD")
D3=ncread(file3, "D")

t4=ncread(file4, "time")/1days
DD4=ncread(file4, "DD")
D4=ncread(file4, "D")

t5=ncread(file5, "time")/1days
DD5=ncread(file5, "DD")
D5=ncread(file5, "D")

t6=ncread(file6, "time")/1days
DD6=ncread(file6, "DD")
D6=ncread(file6, "D")

t7=ncread(file7, "time")/1days
DD7=ncread(file7, "DD")
D7=ncread(file7, "D")

t8=ncread(file8, "time")/1days
DD8=ncread(file8, "DD")
D8=ncread(file8, "D")


V_d = -3.47e-5  #  Detritus sedimentation speed   ms⁻¹
V_dd = -200/day  #  Detritus sedimentation speed  -v_dd_min       50m/day=0.0005878  ms⁻¹ 200/day=2.3e-3m/s
molar_mass_C=12 # 12g/mol
Rd_phy = 6.56 # C:N ratio for P molC molN-1
#G_pump = (D_save[findall(x->x==-98, loc),:]*V_d+DD_save[findall(x->x==-98, loc),:]*V_dd)*Rd_phy*molar_mass_C*1day
z_index = 41 #z1[11]=-112m   new grid z1[41]=-100.43893f0
G_pump0 = (D0[1,1,z_index,:]*V_d+DD0[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump1 = (D1[1,1,z_index,:]*V_d+DD1[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump2 = (D2[1,1,z_index,:]*V_d+DD2[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump3 = (D3[1,1,z_index,:]*V_d+DD3[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump4 = (D4[1,1,z_index,:]*V_d+DD4[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump5 = (D5[1,1,z_index,:]*V_d+DD5[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump6 = (D6[1,1,z_index,:]*V_d+DD6[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump7 = (D7[1,1,z_index,:]*V_d+DD7[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump8 = (D8[1,1,z_index,:]*V_d+DD8[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day

Plots.plot(t0, G_pump0[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump without kelp_base_realforcing",legend=:bottomleft,legendfontsize=5,xlims=(0,730),ylims=(-8000,0),color="gray")
Plots.plot!(t1, G_pump1[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump without kelp_base",legend=:bottomleft,legendfontsize=5,xlims=(0,730),ylims=(-8000,0),color="green")
Plots.plot!(t2, G_pump2[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump kelp density 0.1 continuous",color="blue")
Plots.plot!(t3, G_pump3[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump kelp density 1 continuous",color="red")
Plots.plot!(t4, G_pump4[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump kelp density 10 continuous",color="yellow")
Plots.plot!(t5, G_pump5[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump kelp density 100 continuous",color="black")
Plots.plot!(t6, G_pump6[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump kelp density 1000 pickup",color="orange")
# Plots.plot!(t7, G_pump7[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump kelp density 10000 pickup",color="pink")
# Plots.plot!(t8, G_pump8[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump kelp density 10 pickup",marker = (:circle,1),color="yellow")

# savefig("Gravitational_pump_kelp_com.pdf")


include("old_plotting_uitils.jl")
results0 = load_tracers(file0)
results1 = load_tracers(file1)
results2 = load_tracers(file2)
results3 = load_tracers(file3)
results4 = load_tracers(file4)
results5 = load_tracers(file5)
results6 = load_tracers(file6)
results7 = load_tracers(file7)
results8 = load_tracers(file8)

Nz=50   #33
Lx=Ly=20
pCO₂_0=[]
pCO₂_1=[]
pCO₂_2=[]
pCO₂_3=[]
pCO₂_4=[]
pCO₂_5=[]
pCO₂_6=[]
pCO₂_7=[]
pCO₂_8=[]

t_temperature_node=[0.,65,95,240,364]
temperature_idealize=[9.4,8.1,8.1,13.6,9.1]
temperature_itp = LinearInterpolation((t_temperature_node)days, temperature_idealize) 
t_function(x, y, z, t) = temperature_itp(mod(t, 364days))
s_function(x, y, z, t) = 35.15 

function pco2(result,output)
    for (iter, t) in enumerate(result.t)
        DIC = result.results[7, 1, 1, Nz, iter]   #tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DIC,:ALK)
        ALK = result.results[8, 1, 1, Nz, iter]
        T = t_function(0.5*Lx, 0.5*Ly, 0.0, t) + 273.15
        S = s_function(0.5*Lx, 0.5*Ly, 0.0, t)
        #gas = :CO₂
        #push!(DIC_flux, Boundaries.airseaflux(NaN, NaN, NaN, DIC, ALK, T, S, merge(Boundaries.defaults.airseaflux, (;gas ))))
        push!(output, Boundaries.pCO₂(DIC, ALK, T, S, Boundaries.defaults.airseaflux)) #function pCO₂(DIC, ALK, T, S, params)
    end
end
pco2(results0,pCO₂_0)
pco2(results1,pCO₂_1)
pco2(results2,pCO₂_2)
pco2(results3,pCO₂_3)
pco2(results4,pCO₂_4)
pco2(results5,pCO₂_5)
pco2(results6,pCO₂_6)
pco2(results7,pCO₂_7)
pco2(results8,pCO₂_8)

Plots.plot(t0, pCO₂_0[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "without kelp_base_realforcing",legend=:bottomleft,legendfontsize=5,xlims=(0,730),ylims=(0,500),color="gray")
Plots.plot!(t1, pCO₂_1[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "without kelp_base",legend=:bottomleft,legendfontsize=5,xlims=(0,730),ylims=(0,500),color="green")
Plots.plot!(t2, pCO₂_2[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "kelp density 0.1 continuous",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="blue")
Plots.plot!(t3, pCO₂_3[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "kelp density 1 continuous",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="red")
Plots.plot!(t4, pCO₂_4[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "kelp density 10 continuous",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="yellow")
Plots.plot!(t5, pCO₂_5[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "kelp density 100 continuous",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="black")
# Plots.plot!(t6, pCO₂_6[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "kelp density 1000 pickup",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="orange")
# Plots.plot!(t7, pCO₂_7[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "kelp density 10000 pickup",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="pink")  #pink 
# Plots.plot!(t8, pCO₂_8[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "kelp density 10 pickup",marker = (:circle,1),legend=:bottomleft, xlims=(0,730),ylims=(0,500),color="yellow")   #purple

# pp1 = jldopen(path1*"pco2_water_withoutkelp_long.jld2", "r")
# Plots.plot!(pp1["pCO₂_record"],legend=:bottomleft)
# pp2 = jldopen(path2*"pco2_water_withkelp_density01_largescale_couple.jld2", "r")
# pp3 = jldopen(path3*"pco2_water_withkelp_density1_largescale_couple.jld2", "r")
# pp4 = jldopen(path4*"pco2_water_withkelp_density10_largescale_couple.jld2", "r")
# pp5 = jldopen(path5*"pco2_water_withkelp_density100_largescale_couple.jld2", "r")
# Plots.plot!(366:731,pp2["pCO₂_record"])
# Plots.plot!(366:731,pp3["pCO₂_record"])
# Plots.plot!(366:731,pp4["pCO₂_record"])
# Plots.plot!(366:731,pp5["pCO₂_record"])

savefig("pco2_kelp_com3.pdf")

####################################################################
########### pco2 Mercator 
d1="subpolar_BGC_long.nc"  #Global ocean biogeochemistry hindcast "subpolar_BGC_2010-2020.nc" 2010-1-16-  2020-12-16
# d2="subpolar_BGC_2020-2020forecast.nc"  #Global Ocean Biogeochemistry Analysis and Forecast 2020-1-16-  2020-12-16

spco2_1 = ncread(d1, "spco2");  # in Pa 
spco2_mean = mean(spco2_1, dims=(1,2))[1,1,:]  ## spco2_1_re = reshape(spco2_mean,12,:)  
Plots.plot!(1:730,spco2_mean[1:730]*10,xlabel="Days",ylabel="pCO2 [ppm]",label = "Mercator",marker = (:circle,1),legend=:bottomleft, xlims=(0,730),ylims=(0,500),color="coral") # plot according along column ,starting from Jan 15th ,from pa to ppm 

# savefig("pco2_kelp_com2.pdf")


# #savefig("pco2_data.pdf")

# spco2_2 = ncread(d2, "spco2");  
# spco2_mean2 = mean(spco2_2, dims=(1,2))[1,1,:] 
# # spco2_p2_re = reshape(spco2_mean2) 

# Plots.plot!(1:12,spco2_mean2*10,label = "forecast", xlabel="Months",ylabel="pCO2 [ppm]")
# savefig("pco2_forecast1_hindcast10.pdf")
###########

########### Plot tracers individually   

Rd_phy = 6.56 # C:N ratio for P molC molN-1
phyc1 = ncread(d1, "phyc");  #phyc
phyc1_mean = mean(phyc1, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
depth_phyc1 = ncread(d1, "depth");
P1=Plots.heatmap(1:365, -depth_phyc1[end:-1:1], phyc1_mean[end:-1:1,1:365],interpolate = true,xlabel="Days",ylabel="P [mmolN/m^3]")  #1:730

tracer1=load_tracers(file1)
tracer2=load_tracers(file2)
tracer3=load_tracers(file3)
tracer4=load_tracers(file4)
tracer5=load_tracers(file5)
plot_tracer1=profiles(tracer1)
plot_tracer2=profiles(tracer2)
plot_tracer3=profiles(tracer3)
plot_tracer4=profiles(tracer4)
plot_tracer5=profiles(tracer5)

for (i, tracer) in enumerate(["NO₃", "NH₄", "P", "Z", "D", "DD", "DIC", "ALK", "DOM"])
    Plots.plot(plot_tracer1[i],plot_tracer2[i],plot_tracer3[i],plot_tracer4[i],plot_tracer5[i])
    # savefig(tracer*".pdf")
end

tracer0=load_tracers(file0)
plot_tracer0=profiles(tracer0)
# savefig("realforcing.pdf")


particles5 = load_particles(path5*"particles_density100_largescale_couple.jld2")
plot_particles5 = particles(particles5)
# savefig("kelp_particles_density100.pdf")

# tracer5=load_tracers(file5)
# plot_tracer5=profiles(tracer5)
# title5="kelp_density100.pdf"
# savefig(title5)
