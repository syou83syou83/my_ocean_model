
using NetCDF
using Interpolations
using Statistics
using Plots
using HDF5


filename = "../subpolar.nc" #A small sample of data downloaded is stored in subpolar.nc for ease of use.
# time = ncread(filename, "time")    # time in seconds
# temp = ncread(filename, "temp")    # temperature in Degrees Celsius 
# salinity = ncread(filename, "so")  # salinity in Practical Salinity Unit
# mld = ncread(filename, "mld")      # mixed layer depth in Meters
# par = ncread(filename, "par")      # photosynthetically available radiation in W/m^2

filename1 = "subpolar_physics.nc" # subtropical_physics.nc  subpolar_physics.nc
#ncinfo(filename1)
time_series_second = (0:364)days # start from zero if we don't use extrapolation, we cannot use extrapolation if we wana use annual cycle  
timesave=[time_series_second[i] for i in 1:365]   

so = ncread(filename1, "so");  #salinity
so_scale_factor = ncgetatt(filename1, "so", "scale_factor") #Real_Value = (Display_Value X scale_factor) + add_offset
so_add_offset = ncgetatt(filename1, "so", "add_offset")
salinitysave = mean(so, dims=(1,2))[1:365]*so_scale_factor.+so_add_offset # use [1:365] cause sometimes one year has 366 days. 
# salinity_itp = LinearInterpolation(time_series_second, salinity1) 
#converted to interpolations to access them at arbitary time, how to use it: salinity_itp(mod(timeinseconds,364days))
#plot(salinity)
thetao = ncread(filename1, "thetao");  #temperature
thetao_scale_factor = ncgetatt(filename1, "thetao", "scale_factor") 
thetao_add_offset = ncgetatt(filename1, "thetao", "add_offset")
tempsave = mean(thetao, dims=(1,2))[1:365]*thetao_scale_factor.+thetao_add_offset
# temperature_itp = LinearInterpolation(time_series_second, temperature)  
#plot(temperature)
mlotst = ncread(filename1, "mlotst"); #mixed_layer_depth
mlotst_scale_factor = ncgetatt(filename1, "mlotst", "scale_factor") 
mlotst_add_offset = ncgetatt(filename1, "mlotst", "add_offset")
mldsave = mean(mlotst, dims=(1,2))[1:365]*mlotst_scale_factor.+mlotst_add_offset
# mld_itp = LinearInterpolation(time_series_second, mixed_layer_depth)  


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
parsave=par_mean_timeseries[1,:]


file = "subpolar_Si.nc"
isfile(file) && rm(file)

nccreate(file, "temp", "time", timesave)
ncwrite(tempsave, file, "temp")

nccreate(file, "so", "time", timesave)
ncwrite(salinitysave, file, "so")

nccreate(file, "mld", "time", timesave)
ncwrite(mldsave, file, "mld")

nccreate(file, "par", "time", timesave)
ncwrite(parsave, file, "par")




