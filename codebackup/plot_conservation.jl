
using JLD2
# using Printf
using Oceananigans
using Plots
# using Statistics
using Interpolations
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using NetCDF

include("old_plotting_uitils.jl")

Lx = 20  #20
Ly = 20  #20
Nx = 1
Ny = 1
Nz = 50 # number of points in the vertical direction
Lz = 600 # domain depth

# Generate vertically stretched grid 
refinement = 1.2 # controls spacing near surface (higher means finer spaced)   #10
stretching = 5   # controls rate of stretching at bottom              #5.754
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
file2= path2*"withkelp_density01_largescale_couple.nc"   # 
path3="./continuous_density1_largescale_couple/"
file3= path3*"withkelp_density1_largescale_couple.nc"   # 
path4="./continuous_density10_largescale_couple/"
file4= path4*"withkelp_density10_largescale_couple.nc"   # 
# path5="./continuous_density100_largescale_couple_noflux/"
# file5= path5*"withkelp_density100_largescale_couple_noflux.nc"   #            _noflux/

path5="./temp/"
file5= path5*"withkelp_density100_largescale_couple_noflux.nc"   #  

path6="./pickup_density1000_largescale_couple/"
file6= path6*"withkelp_density1000_largescale.nc"   # 
path7="./pickup_density10000_largescale_couple/"
file7= path7*"withkelp_density10000_largescale.nc"   # 
path8="./pickup_density10_largescale_couple/"
file8= path8*"withkelp_density10_largescale.nc"   # 


############################ plot all tracers in one

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
Z0=ncread(file0, "Z")

NO₃1=ncread(file1, "NO₃")
NH₄1=ncread(file1, "NH₄")
P1=ncread(file1, "P")
D1=ncread(file1, "D")
DD1=ncread(file1, "DD") 
DIC1=ncread(file1, "DIC")
DOM1=ncread(file1, "DOM")
Z1=ncread(file1, "Z")

t2=ncread(file2, "time")/1days  
NO₃2=ncread(file2, "NO₃")
NH₄2=ncread(file2, "NH₄")
P2=ncread(file2, "P")
D2=ncread(file2, "D")
DD2=ncread(file2, "DD")  
DIC2=ncread(file2, "DIC")
DOM2=ncread(file2, "DOM")
Z2=ncread(file2, "Z")

NO₃3=ncread(file3, "NO₃")
NH₄3=ncread(file3, "NH₄")
P3=ncread(file3, "P")
D3=ncread(file3, "D")
DD3=ncread(file3, "DD")  
DIC3=ncread(file3, "DIC")
DOM3=ncread(file3, "DOM")
Z3=ncread(file3, "Z")

NO₃4=ncread(file4, "NO₃")
NH₄4=ncread(file4, "NH₄")
P4=ncread(file4, "P")
D4=ncread(file4, "D")
DD4=ncread(file4, "DD")  
DIC4=ncread(file4, "DIC")
DOM4=ncread(file4, "DOM")
Z4=ncread(file4, "Z")

NO₃5=ncread(file5, "NO₃")
NH₄5=ncread(file5, "NH₄")
P5=ncread(file5, "P")
D5=ncread(file5, "D")
DD5=ncread(file5, "DD")  
DIC5=ncread(file5, "DIC")
DOM5=ncread(file5, "DOM")
Z5=ncread(file5, "Z")

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

######################### check N conservation not done yet 
######################### method 1
volume = grid.Δzᵃᵃᶜ[1:end-3]*grid.Δxᶜᵃᵃ*grid.Δyᵃᶜᵃ
qq0=(NO₃0[1,1,:,:]+NH₄0[1,1,:,:]+P0[1,1,:,:]+D0[1,1,:,:]+DD0[1,1,:,:]+DOM0[1,1,:,:]+Z0[1,1,:,:]).*volume
q0=sum(qq0,dims=1)
qq1=(NO₃1[1,1,:,:]+NH₄1[1,1,:,:]+P1[1,1,:,:]+D1[1,1,:,:]+DD1[1,1,:,:]+DOM1[1,1,:,:]+Z1[1,1,:,:]).*volume
q1=sum(qq1,dims=1)
qq2=(NO₃2[1,1,:,:]+NH₄2[1,1,:,:]+P2[1,1,:,:]+D2[1,1,:,:]+DD2[1,1,:,:]+DOM2[1,1,:,:]+Z2[1,1,:,:]).*volume
q2=sum(qq2,dims=1)
qq3=(NO₃3[1,1,:,:]+NH₄3[1,1,:,:]+P3[1,1,:,:]+D3[1,1,:,:]+DD3[1,1,:,:]+DOM3[1,1,:,:]+Z3[1,1,:,:]).*volume
q3=sum(qq3,dims=1)
qq4=(NO₃4[1,1,:,:]+NH₄4[1,1,:,:]+P4[1,1,:,:]+D4[1,1,:,:]+DD4[1,1,:,:]+DOM4[1,1,:,:]+Z4[1,1,:,:]).*volume
q4=sum(qq4,dims=1)
qq5=(NO₃5[1,1,:,:]+NH₄5[1,1,:,:]+P5[1,1,:,:]+D5[1,1,:,:]+DD5[1,1,:,:]+DOM5[1,1,:,:]+Z5[1,1,:,:]).*volume
q5=sum(qq5,dims=1)
particles2 = load_particles(path2*"particles_density01_largescale_couple.jld2")
particles3 = load_particles(path3*"particles_density1_largescale_couple.jld2")
particles4 = load_particles(path4*"particles_density10_largescale_couple.jld2")
particles5 = load_particles(path5*"particles_density100_largescale_couple_noflux.jld2")
kA=0.6 
density2=0.1
density3=1
density4=10
density5=100

kN=2.72 #g/(gN)^-1

pp2=(0.01.+particles2.results[8,:,:]).*particles2.results[7,:,:]*kA*density2/14*1000   #/kN
p2=sum(pp2,dims=1)
pp3=(0.01.+particles3.results[8,:,:]).*particles3.results[7,:,:]*kA*density3/14*1000
p3=sum(pp3,dims=1)
pp4=(0.01.+particles4.results[8,:,:]).*particles4.results[7,:,:]*kA*density4/14*1000
p4=sum(pp4,dims=1)
pp5=(0.01.+particles5.results[8,:,:]).*particles5.results[7,:,:]*kA*density5/14*1000
p5=sum(pp5,dims=1)

# Plots.plot(q0[:],legend=:topleft,xlabel=xlabel,ylabel="total N (mmolN)",label=title0[1:end-4],title="(NO₃+NH₄+P+D+DD+DOM+Z).*volume")
# Plots.plot!(q1[:],label=title1[1:end-4])
# Plots.plot!(366:731,q2[:]+p2[:],label=title2[1:end-4])
# Plots.plot!(366:731,q3[:]+p3[:],label=title3[1:end-4])
# Plots.plot!(366:731,q4[:]+p4[:],label=title4[1:end-4])
# Plots.plot!(366:731,q5[:]+p5[:],label=title5[1:end-4],color="black")
Plots.plot(q5[:],xlabel=xlabel,ylabel="total N (mmolN)",label="tracers",color="red",legend=:right) #366:731,
Plots.plot!(p5[:],label="kelp",color="blue")
Plots.plot!(q5[:]+p5[:],label="total",color="black")
savefig("N_conservation_80days.pdf")

## Plots.plot!(0:365,q2[:],label=title2[1:end-4])          # check withoutkelp and it's conserved of course it has flux. 
## Plots.plot!(0:365,q3[:],label=title3[1:end-4])
## Plots.plot!(0:365,q5[:],label=title5[1:end-4])
#############################################################



# plot_particles5=particles(particles5)
# savefig("particles5_noflux.pdf")
# julia> particles5.properties
# (:x, :y, :z, :u, :v, :w, :A, :N, :C, :j_NO₃, :j_NH₄, :j_DIC, :j_OXY, :e, :ν, :NO₃, :NH₄, :PAR)




# Plots.plot!(366:731,q5[:]+ppp[:],label=title5[1:end-4])



# tracer5=load_tracers(file5)
# plot_tracer5=profiles(tracer5)
# savefig("kelp_density100_noflux.pdf")





