

using JLD2
# using Printf
# using Oceananigans
using Plots
# using Statistics
using Interpolations
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using NetCDF

filename = "subpolar_Si.nc" #A small sample of data downloaded is stored in subpolar.nc for ease of use.
time = ncread(filename, "time")    # time in seconds
# par = ncread(filename, "par")      # photosynthetically available radiation in W/m^2
# temp = ncread(filename, "temp")    # temperature in Degrees Celsius 
# salinity = ncread(filename, "so")  # salinity in Practical Salinity Unit
mld = ncread(filename, "mld")      # mixed layer depth in Meters

t1_node=[0.,50,80,100,300,364]
mldplus_idealize=[300,500,500,40,40,300] 

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


############################ plot all tracers in one
include("old_plotting_uitils.jl")

tracer0=load_tracers(file0)
plot_tracer0=profiles(tracer0)
title0="realforcing.pdf"
# savefig(title0)
tracer1=load_tracers(file1)
plot_tracer1=profiles(tracer1)
title1="withoutkelpbase_idealize.pdf"
# savefig(title1)
tracer2=load_tracers(file2)
plot_tracer2=profiles(tracer2)
title2="kelp_density01.pdf"
# savefig(title2)
tracer3=load_tracers(file3)
plot_tracer3=profiles(tracer3)
title3="kelp_density1.pdf"
# savefig(title3)
tracer4=load_tracers(file4)
plot_tracer4=profiles(tracer4)
title4="kelp_density10.pdf"
# savefig(title4)
tracer5=load_tracers(file5)
plot_tracer5=profiles(tracer5)
title5="kelp_density100.pdf"
# savefig(title5)


############################ plot one tracers for all cases

z0=ncread(file0, "zC")
t0=ncread(file0, "time")/1days  

NO₃0=ncread(file0, "NO₃")
NH₄0=ncread(file0, "NH₄")
P0=ncread(file0, "P")
D0=ncread(file0, "D")
DD0=ncread(file0, "DD")  
DIC0=ncread(file0, "DIC")
DOM0=ncread(file0, "DOM")

NO₃1=ncread(file1, "NO₃")
NH₄1=ncread(file1, "NH₄")
P1=ncread(file1, "P")
D1=ncread(file1, "D")
DD1=ncread(file1, "DD") 
DIC1=ncread(file1, "DIC")
DOM1=ncread(file1, "DOM")

t2=ncread(file2, "time")/1days  
NO₃2=ncread(file2, "NO₃")
NH₄2=ncread(file2, "NH₄")
P2=ncread(file2, "P")
D2=ncread(file2, "D")
DD2=ncread(file2, "DD")  
DIC2=ncread(file2, "DIC")
DOM2=ncread(file2, "DOM")

NO₃3=ncread(file3, "NO₃")
NH₄3=ncread(file3, "NH₄")
P3=ncread(file3, "P")
D3=ncread(file3, "D")
DD3=ncread(file3, "DD")  
DIC3=ncread(file3, "DIC")
DOM3=ncread(file3, "DOM")

NO₃4=ncread(file4, "NO₃")
NH₄4=ncread(file4, "NH₄")
P4=ncread(file4, "P")
D4=ncread(file4, "D")
DD4=ncread(file4, "DD")  
DIC4=ncread(file4, "DIC")
DOM4=ncread(file4, "DOM")

NO₃5=ncread(file5, "NO₃")
NH₄5=ncread(file5, "NH₄")
P5=ncread(file5, "P")
D5=ncread(file5, "D")
DD5=ncread(file5, "DD")  
DIC5=ncread(file5, "DIC")
DOM5=ncread(file5, "DOM")

fs=4
xlabel="Time (days)"
ylabel="z (m)"
title_logP="log(P)"
title_P="P (mmolN/m^3)"
title_NO₃="NO₃ (mmolN/m^3)"
title_NH₄="NH₄ (mmolN/m^3)"
title_D="D (mmolN/m^3)"
title_DD="DD (mmolN/m^3)"
title_DIC="DIC (mmolC/m^3)"
title_DOM="DOM (mmolN/m^3)"

plot_logP0=Plots.heatmap(t0, z0, log.(P0[1,1,:,:]),xlabel=xlabel,ylabel=ylabel,title=title_logP*title0[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
plot_logP0=Plots.plot!(0:729,-[mld;mld],label = "MLD_Mercator",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,legend=:bottomleft)
plot_logP0=Plots.plot!([t1_node;t1_node.+365],-[mldplus_idealize;mldplus_idealize],label = "MLD_Idealized",legend=:bottomright)

plot_logP1=Plots.heatmap(t0, z0, log.(P1[1,1,:,:]),xlabel=xlabel,ylabel=ylabel,title=title_logP*title1[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
plot_logP1=Plots.plot!(0:729,-[mld;mld],label = "MLD_Mercator",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,legend=:bottomleft)
plot_logP1=Plots.plot!([t1_node;t1_node.+365],-[mldplus_idealize;mldplus_idealize],label = "MLD_Idealized",legend=:bottomright)


plot_logP2=Plots.heatmap(t2, z0, log.(P2[1,1,:,:]),xlabel=xlabel,ylabel=ylabel,title=title_logP*title2[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
plot_logP3=Plots.heatmap(t2, z0, log.(P3[1,1,:,:]),xlabel=xlabel,ylabel=ylabel,title=title_logP*title3[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
plot_logP4=Plots.heatmap(t2, z0, log.(P4[1,1,:,:]),xlabel=xlabel,ylabel=ylabel,title=title_logP*title4[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
plot_logP5=Plots.heatmap(t2, z0, log.(P5[1,1,:,:]),xlabel=xlabel,ylabel=ylabel,title=title_logP*title5[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
Plots.plot(plot_logP0,plot_logP1,plot_logP2,plot_logP3,plot_logP4,plot_logP5)
# savefig("logP.pdf")

plot_P0=Plots.heatmap(t0, z0, P0[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_P*title0[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
plot_P0=Plots.plot!(0:729,-[mld;mld],label = "MLD_Mercator",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,legend=:bottomleft)
plot_P0=Plots.plot!([t1_node;t1_node.+365],-[mldplus_idealize;mldplus_idealize],label = "MLD_Idealized",legend=:bottomright)

plot_P1=Plots.heatmap(t0, z0, P1[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_P*title1[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
plot_P1=Plots.plot!(0:729,-[mld;mld],label = "MLD_Mercator",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,legend=:bottomleft)
plot_P1=Plots.plot!([t1_node;t1_node.+365],-[mldplus_idealize;mldplus_idealize],label = "MLD_Idealized",legend=:bottomright)

plot_P2=Plots.heatmap(t2, z0, P2[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_P*title2[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
plot_P3=Plots.heatmap(t2, z0, P3[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_P*title3[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
plot_P4=Plots.heatmap(t2, z0, P4[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_P*title4[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
plot_P5=Plots.heatmap(t2, z0, P5[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_P*title5[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
Plots.plot(plot_P0,plot_P1,plot_P2,plot_P3,plot_P4,plot_P5)
# savefig("P.pdf")

clims_NO₃=(0,18.5)
plot_NO₃0=Plots.heatmap(t0, z0, NO₃0[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NO₃*title0[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NO₃)  #1:730
plot_NO₃0=Plots.plot!(0:729,-[mld;mld],label = "MLD_Mercator",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,legend=:bottomleft)

plot_NO₃1=Plots.heatmap(t0, z0, NO₃1[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NO₃*title1[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NO₃)  #1:730
plot_NO₃1=Plots.plot!([t1_node;t1_node.+365],-[mldplus_idealize;mldplus_idealize],label = "MLD_Idealized",legend=:bottomright)

plot_NO₃2=Plots.heatmap(t2, z0, NO₃2[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NO₃*title2[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NO₃)  #1:730
plot_NO₃3=Plots.heatmap(t2, z0, NO₃3[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NO₃*title3[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NO₃)  #1:730
plot_NO₃4=Plots.heatmap(t2, z0, NO₃4[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NO₃*title4[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NO₃)  #1:730
plot_NO₃5=Plots.heatmap(t2, z0, NO₃5[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NO₃*title5[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NO₃)  #1:730
Plots.plot(plot_NO₃0,plot_NO₃1,plot_NO₃2,plot_NO₃3,plot_NO₃4,plot_NO₃5)
# savefig("NO₃.pdf")

clims_NH₄=(0,0.5)
plot_NH₄0=Plots.heatmap(t0, z0, NH₄0[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NH₄*title0[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NH₄)  #1:730
plot_NH₄1=Plots.heatmap(t0, z0, NH₄1[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NH₄*title1[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NH₄)  #1:730
plot_NH₄2=Plots.heatmap(t2, z0, NH₄2[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NH₄*title2[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NH₄)  #1:730
plot_NH₄3=Plots.heatmap(t2, z0, NH₄3[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NH₄*title3[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NH₄)  #1:730
plot_NH₄4=Plots.heatmap(t2, z0, NH₄4[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NH₄*title4[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NH₄)  #1:730
plot_NH₄5=Plots.heatmap(t2, z0, NH₄5[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_NH₄*title5[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_NH₄)  #1:730
Plots.plot(plot_NH₄0,plot_NH₄1,plot_NH₄2,plot_NH₄3,plot_NH₄4,plot_NH₄5)
# savefig("NH₄.pdf")


clims_D=(0,0.7)
plot_D0=Plots.heatmap(t0, z0, D0[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_D*title0[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=(0,1.25))  #1:730
plot_D1=Plots.heatmap(t0, z0, D1[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_D*title1[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_D)  #1:730
plot_D2=Plots.heatmap(t2, z0, D2[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_D*title2[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_D)  #1:730
plot_D3=Plots.heatmap(t2, z0, D3[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_D*title3[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_D)  #1:730
plot_D4=Plots.heatmap(t2, z0, D4[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_D*title4[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_D)  #1:730
plot_D5=Plots.heatmap(t2, z0, D5[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_D*title5[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_D)  #1:730
Plots.plot(plot_D0,plot_D1,plot_D2,plot_D3,plot_D4,plot_D5)
# savefig("D.pdf")

clims_DD=(0,0.25)
plot_DD0=Plots.heatmap(t0, z0, DD0[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DD*title0[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=(0,1.7))  #1:730
plot_DD1=Plots.heatmap(t0, z0, DD1[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DD*title1[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DD)  #1:730
plot_DD2=Plots.heatmap(t2, z0, DD2[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DD*title2[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DD)  #1:730
plot_DD3=Plots.heatmap(t2, z0, DD3[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DD*title3[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DD)  #1:730
plot_DD4=Plots.heatmap(t2, z0, DD4[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DD*title4[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DD)  #1:730
plot_DD5=Plots.heatmap(t2, z0, DD5[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DD*title5[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DD)  #1:730
Plots.plot(plot_DD0,plot_DD1,plot_DD2,plot_DD3,plot_DD4,plot_DD5)
# savefig("DD.pdf")

clims_DIC=(1750,2210)
plot_DIC0=Plots.heatmap(t0, z0, DIC0[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DIC*title0[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DIC)  #1:730
plot_DIC1=Plots.heatmap(t0, z0, DIC1[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DIC*title1[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DIC)  #1:730
plot_DIC2=Plots.heatmap(t2, z0, DIC2[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DIC*title2[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DIC)  #1:730
plot_DIC3=Plots.heatmap(t2, z0, DIC3[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DIC*title3[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DIC)  #1:730
plot_DIC4=Plots.heatmap(t2, z0, DIC4[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DIC*title4[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DIC)  #1:730
plot_DIC5=Plots.heatmap(t2, z0, DIC5[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DIC*title5[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DIC)  #1:730
Plots.plot(plot_DIC0,plot_DIC1,plot_DIC2,plot_DIC3,plot_DIC4,plot_DIC5)
# savefig("DIC.pdf")

clims_DOM=(0,2.3)
plot_DOM0=Plots.heatmap(t0, z0, DOM0[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DOM*title0[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DOM)  #1:730
plot_DOM1=Plots.heatmap(t0, z0, DOM1[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DOM*title1[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DOM)  #1:730
plot_DOM2=Plots.heatmap(t2, z0, DOM2[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DOM*title2[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DOM)  #1:730
plot_DOM3=Plots.heatmap(t2, z0, DOM3[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DOM*title3[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DOM)  #1:730
plot_DOM4=Plots.heatmap(t2, z0, DOM4[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DOM*title4[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=clims_DOM)  #1:730
plot_DOM5=Plots.heatmap(t2, z0, DOM5[1,1,:,:],xlabel=xlabel,ylabel=ylabel,title=title_DOM*title5[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=(0,13))  #1:730
Plots.plot(plot_DOM0,plot_DOM1,plot_DOM2,plot_DOM3,plot_DOM4,plot_DOM5)
# savefig("DOM.pdf")


####################################################################
########### plot Mercator 
d1="subpolar_BGC_long.nc"  #Global ocean biogeochemistry hindcast "subpolar_BGC_2010-2020.nc" 2010-1-16-  2020-12-16
# spco2_1 = ncread(d1, "spco2");  # in Pa 
# spco2_mean = mean(spco2_1, dims=(1,2))[1,1,:]  ## spco2_1_re = reshape(spco2_mean,12,:)  
# Plots.plot(spco2_mean[:]*10,label = "Mercator") # plot according along column ,starting from Jan 15th ,from pa to ppm 
#savefig("pco2_data.pdf")

########### Plot tracers individually   

Rd_phy = 6.56 # C:N ratio for P molC molN-1
depth_phyc = ncread(d1, "depth");

phyc = ncread(d1, "phyc");  #phyc
phyc_mean = mean(phyc, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
P_mercator=Plots.heatmap(0:729, -depth_phyc[end:-1:1], phyc_mean[end:-1:1,1:730],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, clims=(0,1.3),xlabel=xlabel,ylabel=ylabel,title="P (mmolN/m^3) Mercator")  #1:730
P_mercator=Plots.plot!(0:729,-[mld;mld],label = "MLD_Mercator",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,legend=:bottomleft)
P_mercator=Plots.plot!([t1_node;t1_node.+365],-[mldplus_idealize;mldplus_idealize],label = "MLD_Idealized",legend=:bottomright)
Plots.plot(P_mercator,plot_P0,plot_P1,layout=(2, 2))
# savefig("P_Mercator.pdf")

logP_mercator=Plots.heatmap(0:729, -depth_phyc[end:-1:1], log.(phyc_mean[end:-1:1,1:730]),titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, xlabel=xlabel,ylabel=ylabel,title="P (mmolN/m^3) Mercator")  #1:730
logP_mercator=Plots.plot!(0:729,-[mld;mld],label = "MLD_Mercator",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,legend=:bottomleft)
logP_mercator=Plots.plot!([t1_node;t1_node.+365],-[mldplus_idealize;mldplus_idealize],label = "MLD_Idealized",legend=:bottomright)
Plots.plot(logP_mercator,plot_logP0,plot_logP1,layout=(2, 2))
# savefig("logP_Mercator.pdf")



no3 = ncread(d1, "no3");  
no3_mean = mean(no3, dims=(1,2))[1,1,:,:]
NO₃_mercator= Plots.heatmap(0:729, -depth_phyc[end:-1:1], no3_mean[end:-1:1,1:730],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, xlabel=xlabel,ylabel=ylabel,title="NO3 (mmolN/m^3) Mercator")
NO₃_mercator=Plots.plot!(0:729,-[mld;mld],label = "MLD_Mercator",titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs,xlabel=xlabel,ylabel=ylabel,legend=:bottomleft)

Plots.plot(NO₃_mercator,plot_NO₃0,plot_NO₃1,layout=(2, 2))
savefig("NO₃_Mercator.pdf")