

using JLD2
# using Printf
# using Oceananigans
using Plots
# using Statistics
using Interpolations
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using NetCDF

file1= "withoutkelp.nc"   #"kelp_example.nc"     
file2= "withkelp.nc"   # "kelp_dense.nc"  
file3= "withoutkelp_long.nc"   # "kelp_dense.nc"  


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

V_d = -3.47e-5  #  Detritus sedimentation speed   ms⁻¹
V_dd = -50/day  #  Detritus sedimentation speed  -v_dd_min       50m/day=0.0005878  ms⁻¹
molar_mass_C=12 # 12g/mol
Rd_phy = 6.56 # C:N ratio for P molC molN-1
#G_pump = (D_save[findall(x->x==-98, loc),:]*V_d+DD_save[findall(x->x==-98, loc),:]*V_dd)*Rd_phy*molar_mass_C*1day
z_index = 13 #z1[11]=-112m
G_pump1 = (D1[1,1,z_index,:]*V_d+DD1[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump2 = (D2[1,1,z_index,:]*V_d+DD2[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
G_pump3 = (D3[1,1,z_index,:]*V_d+DD3[1,1,z_index,:]*V_dd)*Rd_phy*molar_mass_C*1day
Plots.plot(t1, G_pump1[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump without kelp",legend=:bottomright,xlims=(0,730),ylims=(-500,0),color="green")
Plots.plot!(t2, G_pump2[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump with kelp",legend=:bottomright,xlims=(0,730),ylims=(-500,0),color="blue")
Plots.plot!(t3, G_pump3[:], xlabel="Days",ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump without kelp long",legend=:bottomright,xlims=(0,730),ylims=(-500,0),color="red")
savefig("Gravitational_pump_kelp.pdf")


include("old_plotting_uitils.jl")
results1 = load_tracers(file1)
results2 = load_tracers(file2)
results3 = load_tracers(file3)

filename = "./OceanBioME_example_data/subpolar.nc" #A small sample of data downloaded is stored in subpolar.nc for ease of use.
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
s_function(x, y, z, t) = salinity_itp(mod(t, 364days))


Nz=50   #33
Lx=Ly=1
pCO₂_1=[]
pCO₂_2=[]
pCO₂_3=[]
# for (iter, t) in enumerate(results1.t)
#            DIC = results1.results[7, 1, 1, Nz, iter]   #tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DIC,:ALK)
#            ALK = results1.results[8, 1, 1, Nz, iter]
#            T = t_function(0.5*Lx, 0.5*Ly, 0.0, t) + 273.15
#            S = s_function(0.5*Lx, 0.5*Ly, 0.0, t)
#            #gas = :CO₂
#            #push!(DIC_flux, Boundaries.airseaflux(NaN, NaN, NaN, DIC, ALK, T, S, merge(Boundaries.defaults.airseaflux, (;gas ))))
#            push!(pCO₂_1, Boundaries.pCO₂(DIC, ALK, T, S, Boundaries.defaults.airseaflux)) #function pCO₂(DIC, ALK, T, S, params)
# end

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


Plots.plot(t1, pCO₂_1[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "without kelp",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="green")
Plots.plot!(t2, pCO₂_2[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "with kelp",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="blue")
Plots.plot!(t3, pCO₂_3[:], xlabel="Days",ylabel="pCO2 [ppm]",label = "without kelp long",legend=:bottomleft,xlims=(0,730),ylims=(0,500),color="red")

savefig("pco2_kelp.pdf")

########### pco2 from two different database hindcast and forecast
d1="subpolar_BGC_2020-2020hindcast.nc"  #Global ocean biogeochemistry hindcast "subpolar_BGC_2010-2020.nc" 2010-1-16-  2020-12-16
d2="subpolar_BGC_2020-2020forecast.nc"  #Global Ocean Biogeochemistry Analysis and Forecast 2020-1-16-  2020-12-16

spco2_1 = ncread(d1, "spco2");  # in Pa 
spco2_mean = mean(spco2_1, dims=(1,2))[1,1,:]  # 2010-2020 is 11 years monthly data 
spco2_1_re = reshape(spco2_mean,12,:)  

Plots.plot(1:12,spco2_1_re[:,:]*10,label = "hindcast") # plot according along column ,starting from Jan 15th ,from pa to ppm 
#savefig("pco2_data.pdf")

spco2_2 = ncread(d2, "spco2");  
spco2_mean2 = mean(spco2_2, dims=(1,2))[1,1,:] 
# spco2_p2_re = reshape(spco2_mean2) 

Plots.plot!(1:12,spco2_mean2*10,label = "forecast", xlabel="Months",ylabel="pCO2 [ppm]")
savefig("pco2_forecast1_hindcast10.pdf")
###########

########### P from two different database hindcast and forecast

Rd_phy = 6.56 # C:N ratio for P molC molN-1

phyc1 = ncread(d1, "phyc");  #phyc
phyc1_mean = mean(phyc1, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
depth_phyc1 = ncread(d1, "depth");
P1=Plots.heatmap(1:12, -depth_phyc1[end:-1:1], phyc1_mean[end:-1:1,1:12],interpolate = true,xlabel="Months",ylabel="P [mmolN/m^3]")  #1:730

phyc2 = ncread(d2, "phyc");  #phyc
phyc2_mean = mean(phyc2, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
depth_phyc2 = ncread(d2, "depth");
P2=Plots.heatmap(1:12, -depth_phyc2[end:-1:1], phyc2_mean[end:-1:1,1:12],interpolate = true,xlabel="Months",ylabel="P [mmolN/m^3]")  #1:730
P3=Plots.heatmap(1:366, z1, Phyc1[1,1,:,:],interpolate = true,xlabel="Days",ylabel="P [mmolN/m^3]")  #1:730

Plots.plot(P1,P2,P3)
savefig("Heat_P_hind_fore_column.pdf")

p1_itp = LinearInterpolation((z1,1:366), Phyc1[1,1,:,:])  #my simulation
phyc1_itp = LinearInterpolation((depth_phyc1,1:12), phyc1_mean)  #hindcast
phyc2_itp = LinearInterpolation((depth_phyc2,1:12), phyc2_mean)   #forecast

Plots.plot(phyc1_itp(100.0,1:12),label = "hindcast",xlabel="Months",ylabel="P [mmolN/m^3]")
Plots.plot!(phyc2_itp(100.0,1:12),label = "forecast",xlabel="Months",ylabel="P [mmolN/m^3]")

p1_temp=p1_itp(-100.0,1:366)
Plots.plot!((1:366)/30,p1_temp, label = "column",xlabel="Months",ylabel="P [mmolN/m^3]")
savefig("curve_P_100m.pdf")

########### N from two different database hindcast and forecast

n1 = ncread(d1, "no3");  #
n1_mean = mean(n1, dims=(1,2))[1,1,:,:] #from mmolC/m^3 to mmolN/m^3
N1=Plots.heatmap(1:12, -depth_phyc1[end:-1:1], n1_mean[end:-1:1,1:12],interpolate = true,xlabel="Months",ylabel="NO3 [mmolN/m^3]")  #1:730

n2 = ncread(d2, "no3");  #
n2_mean = mean(n2, dims=(1,2))[1,1,:,:] #from mmolC/m^3 to mmolN/m^3
N2=Plots.heatmap(1:12, -depth_phyc2[end:-1:1], n2_mean[end:-1:1,1:12],interpolate = true,xlabel="Months",ylabel="NO3 [mmolN/m^3]")  #1:730
N3=Plots.heatmap(1:366, z1, NO₃1[1,1,:,:],interpolate = true,xlabel="Days",ylabel="NO3 [mmolN/m^3]")  #1:730
Plots.plot(N1,N2,N3)
savefig("Heat_NO3_hind_fore_column.pdf")

N3_itp = LinearInterpolation((z1,1:366), NO₃1[1,1,:,:])  #my simulation
n1_itp = LinearInterpolation((depth_phyc1,1:12), n1_mean)  #hindcast 
n2_itp = LinearInterpolation((depth_phyc2,1:12), n2_mean)   #forecast

Plots.plot(n1_itp(100.0,1:12),label = "hindcast",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
Plots.plot!(n2_itp(100.0,1:12),label = "forecast",xlabel="Months",ylabel="NO3 [mmolN/m^3]")

N3_temp=N3_itp(-100.0,1:366)

Plots.plot!((1:366)/30,N3_temp,label = "column",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
savefig("curve_NO3_100m.pdf")



########### P and N from longer hindcast database
d10="subpolar_BGC_2010-2020.nc"  #Global ocean biogeochemistry hindcast "subpolar_BGC_2010-2020.nc" 2010-1-16-  2020-12-16

phyc10 = ncread(d10, "phyc");  #phyc
phyc10_mean = mean(phyc10, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
phyc10_itp = LinearInterpolation((depth_phyc1,1:132), phyc10_mean)  #hindcast
phyc10_temp=phyc10_itp(100.0,1:132)
phyc10_temp_re=reshape(phyc10_temp,12,:)
Plots.plot(phyc10_temp_re,label = "hindcast",xlabel="Months",ylabel="P [mmolN/m^3]")
savefig("curve_P_longdata.pdf")

phyc10_yrmean=mean(phyc10_temp_re,dims=2)
phyc10_yrstd=std(phyc10_temp_re,dims=2)
Plots.plot(phyc10_yrmean,ribbon=phyc10_yrstd,label = "mean",xlabel="Months",ylabel="P [mmolN/m^3]")
# Plots.plot!(phyc10_yrmean-phyc10_yrstd,label = "mean-std",xlabel="Months",ylabel="P [mmolN/m^3]")
# Plots.plot!(phyc10_yrmean+phyc10_yrstd,label = "mean+std",xlabel="Months",ylabel="P [mmolN/m^3]")
savefig("curve_P_meanstd.pdf")


n10 = ncread(d10, "no3");  #phyc
n10_mean = mean(n10, dims=(1,2))[1,1,:,:] #from mmolC/m^3 to mmolN/m^3
n10_itp = LinearInterpolation((depth_phyc1,1:132), n10_mean)  #hindcast
n10_temp=n10_itp(100.0,1:132)

n10_temp_re=reshape(n10_temp,12,:)
Plots.plot(n10_temp_re,label = "hindcast",xlabel="Months",ylabel="NO3 [mmolN/m^3]",legend=:bottomright)
savefig("curve_NO3_longdata.pdf")


n10_yrmean=mean(n10_temp_re,dims=2)
n10_yrstd=std(n10_temp_re,dims=2)
Plots.plot(n10_yrmean,ribbon=n10_yrstd,label = "mean",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
# Plots.plot!(n10_yrmean-n10_yrstd,label = "mean-std",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
# Plots.plot!(n10_yrmean+n10_yrstd,label = "mean+std",xlabel="Months",ylabel="NO3 [mmolN/m^3]")
savefig("curve_NO3_meanstd.pdf")


Plots.plot(n10_yrmean,ribbon=n10_yrstd,label = "mean",xlabel="Months",ylabel="NO3 [mmolN/m^3]")