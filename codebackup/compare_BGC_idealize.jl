
using JLD2
# using Printf
# using Oceananigans
using Plots
using Statistics
using Interpolations
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using NetCDF


# file1= "withoutkelp.nc"   #"kelp_example.nc"     
# file2= "withkelp.nc"   # "kelp_dense.nc"  
# file3= "withoutkelp_long.nc"   # "kelp_dense.nc"  
# t3=ncread(file3, "time")/1days



#### forecast data 

d1="../subpolar_BGC_long.nc"  #Global ocean biogeochemistry forecast    2020-1-1: 2022-10-7 daily data
spco2_1 = ncread(d1, "spco2");  # in Pa  
spco2_mean = mean(spco2_1, dims=(1,2))[1,1,:] 
spco2_std = std(spco2_1, dims=(1,2))[1,1,:] 
plot(spco2_mean,ribbon=spco2_std,label = "spatial mean",xlabel="Days",ylabel="pCO2 [Pa]",legend=:bottomright)
# savefig("pco2_forecastdata.pdf")

#### simulation 
pco2_sim = jldopen("pco2_water_withoutkelp.jld2", "r")
pco2=pco2_sim["pCO₂_record"]
plot!(1:365,pco2[1:365]/10,label = "column model",right_margin=20Plots.mm)


filename = "../subpolar_Si.nc" #A small sample of data downloaded is stored in subpolar.nc for ease of use.

time = ncread(filename, "time")    # time in seconds
mld = ncread(filename, "mld")      # mixed layer depth in Meters

p1=twinx()
# plot!(pp,time_save/(1day), G_pump[:], ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump",legend=:bottomright,xlims=(0,730),ylims=(-500,0),color="green")

# t1_node=[0.,50,100,130,280,364]
# mldplus_idealize=[350,500,500,60,60,350] 

t1_node=[0.,50,80,100,300,364]
mldplus_idealize=[300,500,500,40,40,300] 

# mld2_itp = LinearInterpolation((t_node)days, mld_idealize)  #in seconds 
plot!(p1,-mld,label = "Mercator MLD",xlabel="Days",ylabel="Mixed Layer Depth (m)",legend=:right,color="red")
plot!(p1,t1_node,-mldplus_idealize,label = "Idealized MLD",xlabel="Days",ylabel="Mixed Layer Depth (m)",legend=:right,color="black")


savefig("pco2_forecastdata_sim3.pdf")



###########

########### P from   forecast and simulation 

Rd_phy = 6.56 # C:N ratio for P molC molN-1

phyc1 = ncread(d1, "phyc");  #phyc
phyc1_mean = mean(phyc1, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
depth_phyc1 = ncread(d1, "depth");
heatmap(1:1011, -depth_phyc1[end:-1:1], phyc1_mean[end:-1:1,1:1011],interpolate = true,xlabel="Days",ylabel="P [mmolN/m^3]")  #1:730

index=23  # about 109.7meter
phyc1_spatialmean = mean(phyc1/Rd_phy, dims=(1,2))[1,1,index,:] #from mmolC/m^3 to mmolN/m^3
phyc1_spatialstd = std(phyc1/Rd_phy, dims=(1,2))[1,1,index,:] #from mmolC/m^3 to mmolN/m^3
plot(phyc1_spatialmean,ribbon=phyc1_spatialstd,label = "spatial mean",xlabel="Days",ylabel="P [mmolN/m^3]",legend=:topright,right_margin=20Plots.mm)
# savefig("P_forecastdata.pdf")



file1= "withoutkelp.nc"   #
z1=ncread(file1, "zC")
index1=40  # z1[40] = -111.24497f0
t1=ncread(file1, "time")/1days
Phyc1=ncread(file1, "P")
NO₃1=ncread(file1, "NO₃")
plot!(1:365,Phyc1[1,1,index1,1:365],label = "column model")

p2=twinx()
plot!(p2,-mld,label = "Mercator MLD",xlabel="Days",ylabel="Mixed Layer Depth (m)",legend=:right,color="red")
plot!(p2,t1_node,-mldplus_idealize,label = "Idealized MLD",xlabel="Days",ylabel="Mixed Layer Depth (m)",legend=:right,color="black")
savefig("P_forecastdata_sim3.pdf")


########### N from  forecast and simulation 

n1 = ncread(d1, "no3");  #
n1_mean = mean(n1, dims=(1,2))[1,1,:,:] #
heatmap(1:1011, -depth_phyc1[end:-1:1], n1_mean[end:-1:1,1:1011],interpolate = true,xlabel="Days",ylabel="NO3 [mmolN/m^3]")  #1:730

n1_spatialmean = mean(n1, dims=(1,2))[1,1,index,:] #from mmolC/m^3 to mmolN/m^3
n1_spatialstd = std(n1, dims=(1,2))[1,1,index,:] #from mmolC/m^3 to mmolN/m^3

plot(n1_spatialmean,ribbon=n1_spatialstd,label = "spatial mean",xlabel="Days",ylabel="NO3 [mmolN/m^3]",legend=:bottomright,right_margin=20Plots.mm)
# savefig("NO3_forecastdata.pdf")

NO₃1=ncread(file1, "NO₃")
plot!(1:365,NO₃1[1,1,index1,1:365],label = "column model")
p3=twinx()
plot!(p3,-mld,label = "Mercator MLD",xlabel="Days",ylabel="Mixed Layer Depth (m)",legend=:right,color="red")
plot!(p3,t1_node,-mldplus_idealize,label = "Idealized MLD",xlabel="Days",ylabel="Mixed Layer Depth (m)",legend=:right,color="black")

savefig("NO3_forecastdata_sim3.pdf")


