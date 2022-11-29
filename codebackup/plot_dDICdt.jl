
using JLD2
# using Printf
# using Oceananigans
using Plots
# using Statistics
using Interpolations
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using NetCDF

# Define the grid
Lx = 20  #20
Ly = 20  #20
Nx = 1
Ny = 1
Nz = 50 # number of points in the vertical direction
Lz = 600 # domain depth

# Generate vertically stretched grid 
refinement = 10 # controls spacing near surface (higher means finer spaced)   #10
stretching = 0.5   # controls rate of stretching at bottom              #5.754
# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz
# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement
# Bottom-intensified stretching function 
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)
grid = RectilinearGrid(size = (Nx, Ny, Nz), 
                       x = (0, Lx),
                       y = (0, Ly),
                       z = z_faces)        #z = z_faces


# path0="./withoutkelpbase_realforcing/"
# file0= path0*"withoutkelp_long_realforcing.nc"   #    
path1="./withoutkelpbase_new/"
# file1= path1*"withoutkelpbase_new.nc"   #    
# path1="./continuous_density10_largescale_couple_new/"
file1DIC= path1*"DIC_withkelp_density10_largescale_couple_new.jld2"   #   
file1pco2= path1*"pco2_water_withkelp_density10_largescale_couple_new.jld2"   #   

path2="./continuous_density10_largescale_couple_new/"
file2DIC= path2*"DIC_withkelp_density10_largescale_couple_new.jld2"   #   
file2pco2= path2*"pco2_water_withkelp_density10_largescale_couple_new.jld2"   #   



# path2="./continuous_density01_largescale_couple/"
# file2= path2*"withkelp_density01_largescale_couple.nc"   # density=0.1
# path3="./continuous_density1_largescale_couple/"
# file3= path3*"withkelp_density1_largescale_couple.nc"   # 
# path4="./continuous_density10_largescale_couple/"
# file4= path4*"withkelp_density10_largescale_couple.nc"   # 
# path5="./continuous_density100_largescale_couple/"
# file5= path5*"withkelp_density100_largescale_couple.nc"   # 
# path6="./pickup_density1000_largescale_couple/"
# file6= path6*"withkelp_density1000_largescale.nc"   # 
# path7="./pickup_density10000_largescale_couple/"
# file7= path7*"withkelp_density10000_largescale.nc"   # 
# path8="./pickup_density10_largescale_couple/"
# file8= path8*"withkelp_density10_largescale.nc"   # 

volume = grid.Δzᵃᵃᶜ[1:end-3]*grid.Δxᶜᵃᵃ*grid.Δyᵃᶜᵃ

dDICdt1 = jldopen(file1DIC, "r")
dDICdt2 = jldopen(file2DIC, "r")
# dDICdt2["z"][8:50]   #0-500meter 
# Plots.plot(sum((dDICdt1["DIC_record"].*volume)[8:50,:],dims=1)[:],ylims=(-10,3),xlabel="Days",ylabel="dDICdt (mmol/s)",label = "dDICdt withoutkelp",legend=:bottomleft,legendfontsize=5,xlims=(0,2200),line=(:solid,:red),right_margin=30Plots.mm)   #(mmol/m^3/s)
sum_dDICdt1 = sum((dDICdt1["DIC_record"].*volume)[8:50,:],dims=1)[:]
sum_dDICdt2 = sum((dDICdt2["DIC_record"].*volume)[8:50,:],dims=1)[:]

Plots.plot(sum_dDICdt1,ylims=(-17,8),xlabel="Days",ylabel="dDICdt (mmol/s)",label = "source+sink withoutkelp",legend=:bottomright,legendfontsize=5,line=(:solid,:red),right_margin=30Plots.mm)   #(mmol/m^3/s)
Plots.plot!(sum_dDICdt2,xlabel="Days",ylabel="dDICdt (mmol/s)",label = "source+sink kelp",legend=:bottomright,legendfontsize=5,line=(:dot,:red),right_margin=30Plots.mm)   #(mmol/m^3/s)

# savefig("dDICdt.pdf")

pco21 = jldopen(file1pco2, "r")
pco22 = jldopen(file2pco2, "r") # pco21["pCO₂_record"]


U_10 = 10.0
pCO2_air = 413.3
flux1 = 7.7e-4*U_10^2*(pco21["pCO₂_record"].-pCO2_air)/years*1000   #(365*24*3600/1000) # from mol/m^2/year to mmol/m^2/s    should be multiply by 1000 to get mmol !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
flux2 = 7.7e-4*U_10^2*(pco22["pCO₂_record"].-pCO2_air)/years*1000   #(365*24*3600/1000) # mmol/m^2/s

area=grid.Δxᶜᵃᵃ*grid.Δyᵃᶜᵃ
# pp=twinx()
sum_flux1 = flux1*area
sum_flux2 = flux2*area
Plots.plot!(sum_dDICdt1+sum_flux1,label = "source+sink+flux withoutkelp",legend=:top,legendfontsize=5,line=(:solid,:blue)) #xlabel="Days",ylabel="CO2 flux (mmol/s)",
Plots.plot!(sum_dDICdt2+sum_flux2,label = "source+sink+flux kelp",legend=:top,legendfontsize=5,line=(:dot,:blue)) #xlabel="Days",ylabel="CO2 flux (mmol/s)",
Plots.plot!([sum_dDICdt1 sum_dDICdt1+sum_flux1 sum_dDICdt2 sum_dDICdt2+sum_flux2],line=[:red :blue :red :blue],linestyle=[:solid :solid :dot :dot],xlims=(450,500),ylims=(-8,-6),tickfontsize=4,inset = (1, bbox(0.1,0.05,0.25,0.25, :bottom, :left)), subplot = 2,legend=false) #label = ["source+sink withoutkelp" "source+sink+flux withoutkelp"],
Plots.plot!([sum_dDICdt1 sum_dDICdt1+sum_flux1 sum_dDICdt2 sum_dDICdt2+sum_flux2],line=[:red :blue :red :blue],linestyle=[:solid :solid :dot :dot],xlims=(1750,1900),ylims=(-0.5,0.25),tickfontsize=4,inset = (1, bbox(0.05,0.75,0.25,0.25, :bottom, :right)), subplot = 3,legend=false) #label = ["source+sink withoutkelp" "source+sink+flux withoutkelp"],

# line=[:red :blue (:dot :red) (:dot :blue)],

# Plots.plot!([sum_dDICdt1+sum_flux1],xlims=(450,500),ylims=(-8,-6),inset = (1, bbox(0.20,0.1,0.4,0.4, :bottom, :left)), subplot = 3,legend=:bottomright,legendfontsize=4)



# Plots.plot!(pp,flux1*area,xlabel="Days",ylabel="CO2 flux (mmol/s)",label = "flux withoutkelp",legend=:bottomright,legendfontsize=5,xlims=(0,2200),line=(:dot,:red))
# Plots.plot!(pp,flux2*area,label = "flux kelp density = 10",line=(:dot,:blue))   
##https://docs.juliaplots.org/latest/gallery/gaston/generated/gaston-ref12/   linestyles 

# savefig("co2_flux.pdf")

# savefig("dDICdt_co2_flux.pdf")


############################ plot all tracers in one
# include("old_plotting_uitils.jl")
# tracer0=load_tracers(file0)
# plot_tracer0=profiles(tracer0)
# title0="realforcing.pdf"
# # savefig(title0)
# tracer1=load_tracers(file1)
# plot_tracer1=profiles(tracer1)
# title1="withoutkelpbase_idealize.pdf"
# # savefig(title1)
# tracer2=load_tracers(file2)
# plot_tracer2=profiles(tracer2)
# title2="kelp_density01.pdf"
# # savefig(title2)
# tracer3=load_tracers(file3)
# plot_tracer3=profiles(tracer3)
# title3="kelp_density1.pdf"
# # savefig(title3)
# tracer4=load_tracers(file4)
# plot_tracer4=profiles(tracer4)
# title4="kelp_density10.pdf"
# # savefig(title4)
# tracer5=load_tracers(file5)
# plot_tracer5=profiles(tracer5)
# title5="kelp_density100.pdf"
# # savefig(title5)
#############################################################

############################ read in all tracers 

# z0=ncread(file0, "zC")
# t0=ncread(file0, "time")/1days  

# NO₃0=ncread(file0, "NO₃")
# NH₄0=ncread(file0, "NH₄")
# P0=ncread(file0, "P")
# D0=ncread(file0, "D")
# DD0=ncread(file0, "DD")  
# DIC0=ncread(file0, "DIC")
# DOM0=ncread(file0, "DOM")

# NO₃1=ncread(file1, "NO₃")
# NH₄1=ncread(file1, "NH₄")
# P1=ncread(file1, "P")
# D1=ncread(file1, "D")
# DD1=ncread(file1, "DD") 
# DIC1=ncread(file1, "DIC")
# DOM1=ncread(file1, "DOM")

# t2=ncread(file2, "time")/1days  
# NO₃2=ncread(file2, "NO₃")
# NH₄2=ncread(file2, "NH₄")
# P2=ncread(file2, "P")
# D2=ncread(file2, "D")
# DD2=ncread(file2, "DD")  
# DIC2=ncread(file2, "DIC")
# DOM2=ncread(file2, "DOM")

# NO₃3=ncread(file3, "NO₃")
# NH₄3=ncread(file3, "NH₄")
# P3=ncread(file3, "P")
# D3=ncread(file3, "D")
# DD3=ncread(file3, "DD")  
# DIC3=ncread(file3, "DIC")
# DOM3=ncread(file3, "DOM")

# NO₃4=ncread(file4, "NO₃")
# NH₄4=ncread(file4, "NH₄")
# P4=ncread(file4, "P")
# D4=ncread(file4, "D")
# DD4=ncread(file4, "DD")  
# DIC4=ncread(file4, "DIC")
# DOM4=ncread(file4, "DOM")

# NO₃5=ncread(file5, "NO₃")
# NH₄5=ncread(file5, "NH₄")
# P5=ncread(file5, "P")
# D5=ncread(file5, "D")
# DD5=ncread(file5, "DD")  
# DIC5=ncread(file5, "DIC")
# DOM5=ncread(file5, "DOM")

# fs=4
# xlabel="Time (days)"
# ylabel="z (m)"
# title_logP="log(P)"
# title_P="P (mmolN/m^3)"
# title_NO₃="NO₃ (mmolN/m^3)"
# title_NH₄="NH₄ (mmolN/m^3)"
# title_D="D (mmolN/m^3)"
# title_DD="DD (mmolN/m^3)"
# title_DIC="DIC (mmolC/m^3)"
# title_DOM="DOM (mmolN/m^3)"

# plot_logP0=Plots.heatmap(t0, z0, log.(P0[1,1,:,:]),xlabel=xlabel,ylabel=ylabel,title=title_logP*title0[1:end-4],titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs)  #1:730
