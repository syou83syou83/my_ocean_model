

using JLD2
# using Printf
# using Oceananigans
using Plots
# using Statistics
using Interpolations
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using NetCDF

# path0="./withkelp/run3/"
# file1= path0*"withoutkelp.nc"   #     
# file2= path0*"withkelp.nc"   # density=0.1

# path1="./continuous_density1_smallscale_couple/"
# file1= path1*"withkelp_density1_smallscale_couple.nc"   #    

path1="./withoutkelpbase/"
file1= path1*"withoutkelp_long.nc"   #    
path2="./continuous_density01_largescale_couple/"
file2= path2*"withkelp_density01_largescale_couple.nc"   # density=0.1
# path3="./continuous_density01_largescale/"
# file3= path3*"withkelp_density01_largescale.nc"   # 
path3="./continuous_density1_largescale_couple/"
file3= path3*"withkelp_density1_largescale_couple.nc"   # 
path4="./continuous_density10_largescale_couple/"
file4= path4*"withkelp_density10_largescale_couple.nc"   # 
path5="./continuous_density100_largescale_couple/"
file5= path5*"withkelp_density100_largescale_couple.nc"   # 
path6="./pickup_density1000_largescale_couple/"
file6= path6*"withkelp_density1000_largescale.nc"   # 
# path8="./continuous_density20_largescale/"
# file8= path8*"withkelp_density20_largescale.nc"   # 
# path9="./pickup_density1_largescale/"
# file9= path9*"withkelp_density1_largescale.nc"   # 
# path10="./pickup_density100_largescale/"
# file10= path10*"withkelp_density100_largescale.nc"   # 



z1=ncread(file1, "zC")
t1=ncread(file1, "time")/1days
DD1=ncread(file1, "DD")  #DD1[1,1,:,:]
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

# t7=ncread(file7, "time")/1days
# DD7=ncread(file7, "DD")
# D7=ncread(file7, "D")

# t8=ncread(file8, "time")/1days
# DD8=ncread(file8, "DD")
# D8=ncread(file8, "D")

# t9=ncread(file9, "time")/1days
# DD9=ncread(file9, "DD")
# D9=ncread(file9, "D")

# t10=ncread(file10, "time")/1days
# DD10=ncread(file10, "DD")
# D10=ncread(file10, "D")

V_d = -3.47e-5  #  Detritus sedimentation speed   ms⁻¹
V_dd = -200/day  #  Detritus sedimentation speed  -v_dd_min       50m/day=0.0005878  ms⁻¹ 200/day=2.3e-3m/s
molar_mass_C=12 # 12g/mol
Rd_phy = 6.56 # C:N ratio for P molC molN-1
#G_pump = (D_save[findall(x->x==-98, loc),:]*V_d+DD_save[findall(x->x==-98, loc),:]*V_dd)*Rd_phy*molar_mass_C*1day
z_index = 41 #z1[11]=-112m   new grid z1[41]=-100.43893f0
G_pump1 = (D1[1,1,z_index,:]*V_d+DD1[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump2 = (D2[1,1,z_index,:]*V_d+DD2[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump3 = (D3[1,1,z_index,:]*V_d+DD3[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump4 = (D4[1,1,z_index,:]*V_d+DD4[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump5 = (D5[1,1,z_index,:]*V_d+DD5[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump6 = (D6[1,1,z_index,:]*V_d+DD6[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
# G_pump7 = (D7[1,1,z_index,:]*V_d+DD7[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
# G_pump8 = (D8[1,1,z_index,:]*V_d+DD8[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
# G_pump9 = (D9[1,1,z_index,:]*V_d+DD9[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
# G_pump10 = (D10[1,1,z_index,:]*V_d+DD10[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day

Plots.plot(t1, G_pump1[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump without kelp_base",legend=:bottomright,xlims=(0,730),ylims=(-1000,0),color="green")
Plots.plot!(t2, G_pump2[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump without kelp first year",legend=:bottomright,xlims=(0,730),ylims=(-1000,0),color="blue")
Plots.plot!(t3, G_pump3[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump continuous__density01_largescale",legend=:bottomright,xlims=(0,730),ylims=(-1000,0),color="red")
Plots.plot!(t4, G_pump4[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump continuous__density1_largescale",legend=:bottomright,xlims=(0,730),ylims=(-1000,0),color="yellow")
Plots.plot!(t5, G_pump5[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump continuous__density10_largescale",legend=:bottomright,xlims=(0,730),ylims=(-1000,0),color="black")
Plots.plot!(t6, G_pump6[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump continuous__density100_largescale",legend=:bottomright,xlims=(0,730),ylims=(-1000,0),color="orange")
# Plots.plot!(t7, G_pump7[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump continuous__density1000_largescale",legend=:bottomright,xlims=(0,730),ylims=(-1000,0),color="pink")
# Plots.plot!(t8, G_pump8[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump continuous_density20_largescale",legend=:bottomright,xlims=(0,730),ylims=(-1000,0),color="purple")
# Plots.plot!(t9, G_pump9[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump pickup__density1_largescale",legend=:bottomright,xlims=(0,730),ylims=(-1000,0),color="yellow", line=(:dot))
# Plots.plot!(t10, G_pump10[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump pickip_density100_largescale",legend=:bottomright,xlims=(0,730),ylims=(-1000,0),color="orange", line=(:dot))

# savefig("Gravitational_pump_kelp_com.pdf")


include("old_plotting_uitils.jl")
results1 = load_tracers(file1)
results2 = load_tracers(file2)
results3 = load_tracers(file3)
results4 = load_tracers(file4)
results5 = load_tracers(file5)
results6 = load_tracers(file6)
# results7 = load_tracers(file7)
# results8 = load_tracers(file8)
# results9 = load_tracers(file9)
# results10 = load_tracers(file10)
# filename = "./OceanBioME_example_data/subpolar.nc" #A small sample of data downloaded is stored in subpolar.nc for ease of use.
filename = "subpolar_Si.nc" #A small sample of data downloaded is stored in subpolar.nc for ease of use.
time = ncread(filename, "time")    # time in seconds
temp = ncread(filename, "temp")    # temperature in Degrees Celsius 
salinity = ncread(filename, "so")  # salinity in Practical Salinity Unit
mld = ncread(filename, "mld")      # mixed layer depth in Meters
par = ncread(filename, "par")      # photosynthetically available radiation in W/m^2

# Linear interpolation to access temperature, salinity, mld, and surface PAR at arbitrary time
temperature_itp = LinearInterpolation(time, temp) 
salinity_itp = LinearInterpolation(time, salinity) 
mld_itp = LinearInterpolation(time, mld) 
PAR_itp = LinearInterpolation(time, par)

# Define temperature and salinity as functions of x, y, z, and t(in seconds). The temperature and salinity functions are needed to calculate the air-sea CO2 flux.
t_function(x, y, z, t) = temperature_itp(mod(t, 364days)) # the remainder of t after floored division by 364days. It creates an annual cycle representation of temperature. 
s_function(x, y, z, t) = 35.15 # salinity_itp(mod(t, 364days))


Nz=50   #33
Lx=Ly=1
pCO₂_1=[]
pCO₂_2=[]
pCO₂_3=[]
pCO₂_4=[]
pCO₂_5=[]
pCO₂_6=[]
pCO₂_7=[]
pCO₂_8=[]
pCO₂_9=[]
pCO₂_10=[]
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


pco2(results1,pCO₂_1)
pco2(results2,pCO₂_2)
pco2(results3,pCO₂_3)
pco2(results4,pCO₂_4)
pco2(results5,pCO₂_5)
pco2(results6,pCO₂_6)
# pco2(results7,pCO₂_7)
# pco2(results8,pCO₂_8)
# pco2(results9,pCO₂_9)
# pco2(results10,pCO₂_10)

Plots.plot(t1, pCO₂_1[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "without kelp_base",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="green")
Plots.plot!(t2, pCO₂_2[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "with kelp first year",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="blue")
Plots.plot!(t3, pCO₂_3[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "with kelp continuous__density01_largescale",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="red")
Plots.plot!(t4, pCO₂_4[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "with kelp continuous__density1_largescale",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="yellow")
Plots.plot!(t5, pCO₂_5[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "with kelp continuous__density10_largescale",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="black")
Plots.plot!(t6, pCO₂_6[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "with kelp continuous__density100_largescale",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="orange")
# Plots.plot!(t7, pCO₂_7[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "with kelp continuous__density1000_largescale",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="pink")
# Plots.plot!(t8, pCO₂_8[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "with kelp continuous_density20_largescale",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="purple")
# Plots.plot!(t9, pCO₂_9[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "with kelp pickup_density1_largescale",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="yellow", line=(:dot))
# Plots.plot!(t10, pCO₂_10[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "with kelp pickup_density100_largescale",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="orange",line=(:dot))

# Plots.plot!(t3, pCO₂_3[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "without kelp long",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="red")

# savefig("pco2_kelp_com.pdf")


# pp1 = jldopen(path1*"pco2_water_withkelp_density1_smallscale_couple.jld2", "r")
# pp1["pCO₂_record"]
# pp10 = jldopen(path2*"pco2_water_withkelp_density10.jld2", "r")
# pp10["pCO₂_record"]


####################################################################
########### pco2 from two different database hindcast and forecast
d1="subpolar_BGC_long.nc"  #Global ocean biogeochemistry hindcast "subpolar_BGC_2010-2020.nc" 2010-1-16-  2020-12-16
# d2="subpolar_BGC_2020-2020forecast.nc"  #Global Ocean Biogeochemistry Analysis and Forecast 2020-1-16-  2020-12-16

# spco2_1 = ncread(d1, "spco2");  # in Pa 
# spco2_mean = mean(spco2_1, dims=(1,2))[1,1,:]  # 2010-2020 is 11 years monthly data 
# spco2_1_re = reshape(spco2_mean,12,:)  

# Plots.plot(1:12,spco2_1_re[:,:]*10,label = "hindcast") # plot according along column ,starting from Jan 15th ,from pa to ppm 
# #savefig("pco2_data.pdf")

# spco2_2 = ncread(d2, "spco2");  
# spco2_mean2 = mean(spco2_2, dims=(1,2))[1,1,:] 
# # spco2_p2_re = reshape(spco2_mean2) 

# Plots.plot!(1:12,spco2_mean2*10,label = "forecast", xlabel="Months",ylabel="pCO2 [ppm]")
# savefig("pco2_forecast1_hindcast10.pdf")
###########

########### P  

Rd_phy = 6.56 # C:N ratio for P molC molN-1

phyc1 = ncread(d1, "phyc");  #phyc
phyc1_mean = mean(phyc1, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
depth_phyc1 = ncread(d1, "depth");
# P1=Plots.heatmap(1:365, -depth_phyc1[end:-1:1], phyc1_mean[end:-1:1,1:365],interpolate = true,xlabel="Days",ylabel="P [mmolN/m^3]")  #1:730

# phyc2 = ncread(d2, "phyc");  #phyc
# phyc2_mean = mean(phyc2, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
# depth_phyc2 = ncread(d2, "depth");
# P2=Plots.heatmap(1:12, -depth_phyc2[end:-1:1], phyc2_mean[end:-1:1,1:12],interpolate = true,xlabel="Months",ylabel="P [mmolN/m^3]")  #1:730
# P3=Plots.heatmap(1:365, z1, Phyc1[1,1,:,1:365],interpolate = true,xlabel="Days",ylabel="P [mmolN/m^3]")  #1:730

# Plots.plot(P1,P3)
# savefig("P_Mercator_Sim.pdf")


# p1=load_particles(path3*"particles_density01_largescale.jld2")
# plot_part01=particles(p1)
# savefig("p_01.pdf")

# t1=load_tracers(path3*"withkelp_density01_largescale.nc")
# plot_tracer01=profiles(t1)
# savefig("t_01.pdf")

# p6=load_particles(path6*"particles_density100_largescale.jld2")
# plot_part6=particles(p6)
# savefig("p_100.pdf")

# t1=load_tracers(path1*"withoutkelp_density1_smallscale_couple.nc")
# plot_tracer1=profiles(t1)
# savefig("withoutkelp_1.pdf")





# p1_itp = LinearInterpolation((z1,1:366), Phyc1[1,1,:,:])  #my simulation
# phyc1_itp = LinearInterpolation((depth_phyc1,1:12), phyc1_mean)  #hindcast
# phyc2_itp = LinearInterpolation((depth_phyc2,1:12), phyc2_mean)   #forecast

# Plots.plot(phyc1_itp(100.0,1:12),label = "hindcast",xlabel="Months",ylabel="P [mmolN/m^3]")
# Plots.plot!(phyc2_itp(100.0,1:12),label = "forecast",xlabel="Months",ylabel="P [mmolN/m^3]")

# p1_temp=p1_itp(-100.0,1:366)
# Plots.plot!((1:366)/30,p1_temp, label = "column",xlabel="Months",ylabel="P [mmolN/m^3]")
# savefig("curve_P_100m.pdf")

# ########### N from two different database hindcast and forecast

# n1 = ncread(d1, "no3");  #
# n1_mean = mean(n1, dims=(1,2))[1,1,:,:] #from mmolC/m^3 to mmolN/m^3
# N1=Plots.heatmap(1:12, -depth_phyc1[end:-1:1], n1_mean[end:-1:1,1:12],interpolate = true,xlabel="Months",ylabel="NO3 [mmolN/m^3]")  #1:730

# n2 = ncread(d2, "no3");  #
# n2_mean = mean(n2, dims=(1,2))[1,1,:,:] #from mmolC/m^3 to mmolN/m^3
# N2=Plots.heatmap(1:12, -depth_phyc2[end:-1:1], n2_mean[end:-1:1,1:12],interpolate = true,xlabel="Months",ylabel="NO3 [mmolN/m^3]")  #1:730
# N3=Plots.heatmap(1:366, z1, NO₃1[1,1,:,:],interpolate = true,xlabel="Days",ylabel="NO3 [mmolN/m^3]")  #1:730
# Plots.plot(N1,N2,N3)
# savefig("Heat_NO3_hind_fore_column.pdf")

# N3_itp = LinearInterpolation((z1,1:366), NO₃1[1,1,:,:])  #my simulation
# n1_itp = LinearInterpolation((depth_phyc1,1:12), n1_mean)  #hindcast 
# n2_itp = LinearInterpolation((depth_phyc2,1:12), n2_mean)   #forecast

# Plots.plot(n1_itp(100.0,1:12),label = "hindcast",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
# Plots.plot!(n2_itp(100.0,1:12),label = "forecast",xlabel="Months",ylabel="NO3 [mmolN/m^3]")

# N3_temp=N3_itp(-100.0,1:366)

# Plots.plot!((1:366)/30,N3_temp,label = "column",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
# savefig("curve_NO3_100m.pdf")



# ########### P and N from longer hindcast database
# d10="subpolar_BGC_2010-2020.nc"  #Global ocean biogeochemistry hindcast "subpolar_BGC_2010-2020.nc" 2010-1-16-  2020-12-16

# phyc10 = ncread(d10, "phyc");  #phyc
# phyc10_mean = mean(phyc10, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
# phyc10_itp = LinearInterpolation((depth_phyc1,1:132), phyc10_mean)  #hindcast
# phyc10_temp=phyc10_itp(100.0,1:132)
# phyc10_temp_re=reshape(phyc10_temp,12,:)
# Plots.plot(phyc10_temp_re,label = "hindcast",xlabel="Months",ylabel="P [mmolN/m^3]")
# savefig("curve_P_longdata.pdf")

# phyc10_yrmean=mean(phyc10_temp_re,dims=2)
# phyc10_yrstd=std(phyc10_temp_re,dims=2)
# Plots.plot(phyc10_yrmean,ribbon=phyc10_yrstd,label = "mean",xlabel="Months",ylabel="P [mmolN/m^3]")
# # Plots.plot!(phyc10_yrmean-phyc10_yrstd,label = "mean-std",xlabel="Months",ylabel="P [mmolN/m^3]")
# # Plots.plot!(phyc10_yrmean+phyc10_yrstd,label = "mean+std",xlabel="Months",ylabel="P [mmolN/m^3]")
# savefig("curve_P_meanstd.pdf")


# n10 = ncread(d10, "no3");  #phyc
# n10_mean = mean(n10, dims=(1,2))[1,1,:,:] #from mmolC/m^3 to mmolN/m^3
# n10_itp = LinearInterpolation((depth_phyc1,1:132), n10_mean)  #hindcast
# n10_temp=n10_itp(100.0,1:132)

# n10_temp_re=reshape(n10_temp,12,:)
# Plots.plot(n10_temp_re,label = "hindcast",xlabel="Months",ylabel="NO3 [mmolN/m^3]",legend=:bottomright)
# savefig("curve_NO3_longdata.pdf")


# n10_yrmean=mean(n10_temp_re,dims=2)
# n10_yrstd=std(n10_temp_re,dims=2)
# Plots.plot(n10_yrmean,ribbon=n10_yrstd,label = "mean",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
# # Plots.plot!(n10_yrmean-n10_yrstd,label = "mean-std",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
# # Plots.plot!(n10_yrmean+n10_yrstd,label = "mean+std",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
# savefig("curve_NO3_meanstd.pdf")


# Plots.plot(n10_yrmean,ribbon=n10_yrstd,label = "mean",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
# ##################################################################


