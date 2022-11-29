
using NetCDF
using Interpolations
using Statistics
using Plots
using HDF5


using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

# filename = "../subpolar.nc" #A small sample of data downloaded is stored in subpolar.nc for ease of use.
filename = "subpolar_Si.nc" #A small sample of data downloaded is stored in subpolar.nc for ease of use.

time = ncread(filename, "time")    # time in seconds
temp = ncread(filename, "temp")    # temperature in Degrees Celsius 
salinity = ncread(filename, "so")  # salinity in Practical Salinity Unit
mld = ncread(filename, "mld")      # mixed layer depth in Meters
par = ncread(filename, "par")      # photosynthetically available radiation in W/m^2

filename1 = "subpolar_physics.nc" # subtropical_physics.nc  subpolar_physics.nc
#ncinfo(filename1)
time_series_second = (0:364)days # start from zero if we don't use extrapolation, we cannot use extrapolation if we wana use annual cycle  
# time=[time_series_second[i] for i in 1:365]   

so = ncread(filename1, "so");  #salinity
so_scale_factor = ncgetatt(filename1, "so", "scale_factor") #Real_Value = (Display_Value X scale_factor) + add_offset
so_add_offset = ncgetatt(filename1, "so", "add_offset")
salinity1 = mean(so, dims=(1,2))[1:365]*so_scale_factor.+so_add_offset # use [1:365] cause sometimes one year has 366 days. 
salinity_itp = LinearInterpolation(time_series_second, salinity1) 
#converted to interpolations to access them at arbitary time, how to use it: salinity_itp(mod(timeinseconds,364days))
#plot(salinity)
thetao = ncread(filename1, "thetao");  #temperature
thetao_scale_factor = ncgetatt(filename1, "thetao", "scale_factor") 
thetao_add_offset = ncgetatt(filename1, "thetao", "add_offset")
temperature = mean(thetao, dims=(1,2))[1:365]*thetao_scale_factor.+thetao_add_offset
temperature_itp = LinearInterpolation(time_series_second, temperature)  
#plot(temperature)
mlotst = ncread(filename1, "mlotst"); #mixed_layer_depth
mlotst_scale_factor = ncgetatt(filename1, "mlotst", "scale_factor") 
mlotst_add_offset = ncgetatt(filename1, "mlotst", "add_offset")
mixed_layer_depth = mean(mlotst, dims=(1,2))[1:365]*mlotst_scale_factor.+mlotst_add_offset
mld_itp = LinearInterpolation(time_series_second, mixed_layer_depth)  


path="../../examples/OceanBioME_example_data/subpolar/"    #subtropical   #./subpolar/
par_mean_timeseries=zeros(1,365)
for i in 1:365    #https://discourse.julialang.org/t/leading-zeros/30450
    string_i = lpad(string(i), 3, '0')
    filename3=path*"V2020"*string_i*".L3b_DAY_SNPP_PAR.x.nc"
    fid = h5open(filename3, "r")
    par1=read(fid["level-3_binned_data/par"])
    BinList=read(fid["level-3_binned_data/BinList"])  #(:bin_num, :nobs, :nscenes, :weights, :time_rec) 
    par_mean_timeseries[1,i] = mean([par1[i][1]/BinList[i][4] for i in 1:length(par1)])*3.99e-10*545e12/(1day)  #from einstin/m^2/day to W/m^2
end
#a=par_mean_timeseries[1,:]





########## design idealized piecewiselinear interplation mld function 
t_node=[0.,50,86,105,280,364] # in days specify piecewiselinear nodes, better include boundary nodes[0.0,364], must be float 
mld_idealize=[250,420,420,20,20,250] # can either extract from Mercator model mld_itp.(t) or specify positive mld at above nodes

t1_node=[0.,50,100,130,280,364]
mldplus_idealize=[350,500,500,60,60,350] 
mld2_itp = LinearInterpolation((t_node)days, mld_idealize)  #in seconds 
plot(-mld)
plot!(t_node,-mld_idealize)
savefig("mld_idealize.pdf")

t_temperature_node=[0.,65,95,240,364]
temperature_idealize=[9.4,8.1,8.1,13.6,9.1]
temperature_itp = LinearInterpolation((t_temperature_node)days, temperature_idealize) 
plot(temp)
plot!(t_temperature_node,temperature_idealize)
# savefig=("temp_idealize.pdf")

t_par_node=[0.,30,120,200,330,364]
par_idealize=[3,3,90,90,3,3]
plot(par)
plot!(t_par_node,par_idealize)
# savefig=("par_idealize.pdf")



# filename_database = "subpolar_BGC.nc"   #subpolar_BGC
# Rd_phy = 6.56 # C:N ratio for P molC molN-1
# phyc = ncread(filename_database, "phyc");  #phyc
# phyc_mean = mean(phyc, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
# depth_phyc = ncread(filename_database, "depth");
# heatmap(0:729, -depth_phyc[end:-1:1], log.(phyc_mean[end:-1:1,1:730]),interpolate = true)
# plot!([t_node;t_node.+365], -[mld_idealize;mld_idealize],xlabel="time (days)", ylabel="z (m)",xlims=(0,730),ylims=(-600,0),color="blue",right_margin=20Plots.mm)
# plot!([t1_node;t1_node.+365], -[mldplus_idealize;mldplus_idealize],xlabel="time (days)", ylabel="z (m)",xlims=(0,730),ylims=(-600,0),color="green",right_margin=20Plots.mm)

# κ(z, t) = 1e-2*max(1-(z+mld2_itp(mod(t,364days))/2)^2/(mld2_itp(mod(t,364days))/2)^2,0)+1e-4; #setup viscosity and diffusivity in the following Model instantiation
