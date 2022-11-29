#####################in postprocessing to calculate dDICdt instead of using callback-obtained dDICdt 
using JLD2
# using Printf
# using Oceananigans
using Plots
# using Statistics
using Interpolations
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using NetCDF

include("old_plotting_uitils.jl")
duration = 6years

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

t_temperature_node=[0.,66,95,240,364]
temperature_idealize=[9.0,8.05,8.05,13.65,9.0]
temperature_itp = LinearInterpolation((t_temperature_node)days, temperature_idealize) 
t_function(x, y, z, t) = temperature_itp(mod(t, 364days))
s_function(x, y, z, t) = 35.15 

# path0="./withoutkelpbase_realforcing/"
# file0= path0*"withoutkelp_long_realforcing.nc"   #    
path1="./withoutkelpbase_new/"                             # base 
file1= path1*"withoutkelpbase_new.nc"        
 
path2="./continuous_density10_largescale_couple_new/"             #kelp density = 10
file2= path2*"withkelp_density10_largescale_couple_new.nc"   #  

path3="./withoutkelpbase_new3_2/"                    #DIC=2226 1-6 yr
file3= path3*"withoutkelpbase_new3_2.nc"   # 
path4="./withoutkelpbase_new3_2_pickup/"             #DIC=2226 6-12 yr
file4= path4*"withoutkelpbase_new3_2_pickup.nc"   # 

path5="./withoutkelpbase_new2_4/"                    #pco2air=400 1-6 yr
file5= path5*"withoutkelpbase_new2_4.nc"   # 

path6="./withoutkelpbase_new2_pickup/"                    #base 6-12yrs
file6= path6*"withoutkelpbase_new2_pickup.nc"   # 

path7="./continuous_density10_largescale_couple_new_pickup/"             #kelp density = 10   ,6-12yrs
file7= path7*"withkelp_density10_largescale_couple_new_pickup.nc"   #  
path8="./withoutkelpbase_new2_4_pickup/"                    #pco2air=400 1-6 yr
file8= path8*"withoutkelpbase_new2_4_pickup.nc"   # 
# path3="./test/"
# file4DIC= path4*"DIC_withoutkelpbase_new3_2_pickup.jld2"   #   
# file4flux= path4*"fluxco2_withoutkelpbase_new3_2_pickup.jld2"   #   
# file4pco2= path4*"pco2_water_withoutkelpbase_new3_2_pickup.jld2"   #   
# DIC44 = jldopen(file4DIC, "r")["DIC_record"]
# flux44 = jldopen(file4flux, "r")["fluxco2_record"]
# pco244 = jldopen(file4pco2, "r")["pCO₂_record"]

# diffusion = jldopen(file3diff, "r")["diffusion_record"]
# diffusion1 = jldopen(file3diff1, "r")["diffusion1_record"]
# diffusion2 = jldopen(file3diff2, "r")["diffusion2_record"]
# diffusion3 = jldopen(file3diff3, "r")["diffusion3_record"]


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
area=grid.Δxᶜᵃᵃ*grid.Δyᵃᶜᵃ

function moving_avg(arr)
    return [sum(arr[i:i+364]) for i in 1:length(arr)-365]*1day/365days
end
function moving_avg_percent(arr)  
    tem=moving_avg(arr)
    return [tem[i+365]/tem[i]-1 for i in 1:length(tem)-365]
end


results1 = load_tracers(file1)
results2 = load_tracers(file2)
results3 = load_tracers(file3)
results4 = load_tracers(file4)
results5 = load_tracers(file5)
results6 = load_tracers(file6)
results7 = load_tracers(file7)
results8 = load_tracers(file8)

#######################################calculate co2 flux 
function flux(result,output)
    for (iter, t) in enumerate(result.t)
        DIC = result.results[7, 1, 1, Nz, iter]   #tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DIC,:ALK)
        ALK = result.results[8, 1, 1, Nz, iter]
        T = t_function(0.5*Lx, 0.5*Ly, 0.0, t)   # here there is no 273.15 for flux , but pco2 calculation has 273.15
        S = s_function(0.5*Lx, 0.5*Ly, 0.0, t)
        gas = :CO₂
        push!(output, Boundaries.airseaflux(NaN, NaN, NaN, DIC, ALK, T, S, merge(Boundaries.defaults.airseaflux, (;gas ))))
        # push!(output, Boundaries.airseaflux(0.5*Lx,0.5*Ly, t, DIC, ALK, T, S, Boundaries.defaults.airseaflux)) #function pCO₂(DIC, ALK, T, S, params)
    end
end
flux_1=[]
flux_2=[]
flux_3=[]
flux_4=[]
flux_5=[]
flux_6=[]
flux_7=[]
flux_8=[]

flux(results1,flux_1)
flux(results2,flux_2)
flux(results3,flux_3)
flux(results4,flux_4)
flux(results5,flux_5)
flux(results6,flux_6)
flux(results7,flux_7)
flux(results8,flux_8)

##############################################
############################################## calculate pco2 

function pco2(result,output)
    for (iter, t) in enumerate(result.t)
        DIC = result.results[7, 1, 1, Nz, iter]   #tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DIC,:ALK)
        ALK = result.results[8, 1, 1, Nz, iter]
        T = t_function(0.5*Lx, 0.5*Ly, 0.0, t) + 273.15   # here there is no 273.15 for flux , but has 273.15 for pco2 calculation 
        S = s_function(0.5*Lx, 0.5*Ly, 0.0, t)
        push!(output, Boundaries.pCO₂(DIC, ALK, T, S, Boundaries.defaults.airseaflux))
        # push!(output, Boundaries.airseaflux(0.5*Lx,0.5*Ly, t, DIC, ALK, T, S, Boundaries.defaults.airseaflux)) #function pCO₂(DIC, ALK, T, S, params)
    end
end
pco2_1=[]
pco2_2=[]
pco2_3=[]
pco2_4=[]
pco2_5=[]
pco2_6=[]
pco2_7=[]
pco2_8=[]

pco2(results1,pco2_1)
pco2(results2,pco2_2)
pco2(results3,pco2_3)
pco2(results4,pco2_4)
pco2(results5,pco2_5)
pco2(results6,pco2_6)
pco2(results7,pco2_7)
pco2(results8,pco2_8)

##############################################
############################################## calculate DIC forcing   
function DICforcing(result,output)
    for (iter, t) in enumerate(result.t)
        #tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DIC, :ALK, :DOM, :OXY,:PAR)
        for j in 1:Nz
            NO₃ = result.results[1, 1, 1, j, iter] 
            NH₄ = result.results[2, 1, 1, j, iter] 
            P = result.results[3, 1, 1, j, iter] 
            Z = result.results[4, 1, 1, j, iter] 
            D = result.results[5, 1, 1, j, iter] 
            DD = result.results[6, 1, 1, j, iter] 
            DIC = result.results[7, 1, 1, j, iter]   
            ALK = result.results[8, 1, 1, j, iter]
            DOM = result.results[9, 1, 1, j, iter]
            PAR = result.results[11, 1, 1, j, iter]
            # output[j,floor(Int,t/days)+1] = LOBSTER.DIC_forcing(0.5*Lx, 0.5*Ly, result.z[j], t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR, LOBSTER.defaults)
            output[j,iter] = LOBSTER.DIC_forcing(0.5*Lx, 0.5*Ly, result.z[j], t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR, LOBSTER.defaults)

            # DIC_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params)
        end
    end
end

DICforcing_1 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results1,DICforcing_1)
DICforcing_2 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results2,DICforcing_2)
DICforcing_3 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results3,DICforcing_3)
DICforcing_4 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results4,DICforcing_4)
DICforcing_5 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results5,DICforcing_5)
DICforcing_6 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results6,DICforcing_6)
DICforcing_7 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results7,DICforcing_7)
DICforcing_8 = zeros(Nz,floor(Int,duration/days)+1)
DICforcing(results8,DICforcing_8)
##############################################
############################################## plot

fluxsum1 = flux_1*area
fluxsum2 = flux_2*area
fluxsum3 = flux_3*area
fluxsum4 = flux_4*area
fluxsum5 = flux_5*area
fluxsum6 = flux_6*area
fluxsum7 = flux_7*area
fluxsum8 = flux_8*area

fluxmean1 = [mean(fluxsum1[(1:365).+x*365]) for x in [0:5;]]
fluxmean2 = [mean(fluxsum2[(1:365).+x*365]) for x in [0:5;]]
fluxmean3 = [mean(fluxsum3[(1:365).+x*365]) for x in [0:5;]]
fluxmean4 = [mean(fluxsum4[(1:365).+x*365]) for x in [0:5;]]
fluxmean5 = [mean(fluxsum5[(1:365).+x*365]) for x in [0:5;]]
fluxmean6 = [mean(fluxsum6[(1:365).+x*365]) for x in [0:5;]]
fluxmean7 = [mean(fluxsum7[(1:365).+x*365]) for x in [0:4;]]
fluxmean8 = [mean(fluxsum8[(1:365).+x*365]) for x in [0:5;]]


depth_index=1:Nz   #1:Nz 

DICforcingsum1 = sum((DICforcing_1.*volume)[depth_index,:],dims=1)[:]  
DICforcingsum2 = sum((DICforcing_2.*volume)[depth_index,:],dims=1)[:]  
DICforcingsum3 = sum((DICforcing_3.*volume)[depth_index,:],dims=1)[:]  
DICforcingsum4 = sum((DICforcing_4.*volume)[depth_index,:],dims=1)[:]  
DICforcingsum5 = sum((DICforcing_5.*volume)[depth_index,:],dims=1)[:]  
DICforcingsum6 = sum((DICforcing_6.*volume)[depth_index,:],dims=1)[:]  
DICforcingsum7 = sum((DICforcing_7.*volume)[depth_index,:],dims=1)[:]  
DICforcingsum8 = sum((DICforcing_8.*volume)[depth_index,:],dims=1)[:]  

## Plots.plot(DICforcingsum1,ylabel="(mmol/s)",line=(:dot,:blue),ylims=(-0.6,0.3),label="integrated DIC forcing without kelp",legend=:topleft, right_margin=30Plots.mm)
## Plots.plot(DICforcingsum2,ylabel="(mmol/s)",line=(:dot,:blue),ylims=(-0.5,0.3),label="integrated DIC forcing kelp",legend=:topleft, right_margin=30Plots.mm)

# Plots.plot(DICforcingsum1-fluxsum1,line=(:solid,:black),ylabel="(mmol/s)",ylims=(-0.4,0.5),label="integrated DIC forcing+flux without kelp",legend=:topleft,right_margin=30Plots.mm)
# Plots.plot!(DICforcingsum2-fluxsum2,line=(:dot,:black),label="integrated DIC forcing+flux kelp")

# Plots.plot(DICforcingsum1,line=(:solid,:black),ylabel="(mmol/s)",ylims=(-0.5,0.4),label="integrated DIC forcing without kelp",legend=:topleft)
# Plots.plot!(DICforcingsum2,line=(:dot,:black),ylabel="(mmol/s)",ylims=(-0.5,0.4),legend=:topleft,label="integrated DIC forcing kelp")

Plots.plot([results1.t/days;results6.t/days],-[fluxsum1;fluxsum6],line=(:dot,:red),ylabel="(mmol/s)",ylims=(-0.1,0.2),legend=:bottomright,legendfontsize=5,label="integrated flux without kelp") #, right_margin=30Plots.mm
Plots.plot!([results2.t/days;results7.t/days],-[fluxsum2;fluxsum7],line=(:dot,:black),label="integrated flux kelp")
Plots.plot!([results3.t/days;results4.t/days],-[fluxsum3;fluxsum4],line=(:dot,:green),label="integrated flux withoutkelp DIC=2226")
Plots.plot!([results5.t/days;results8.t/days],-[fluxsum5;fluxsum8],line=(:dot,:blue),label="integrated flux withoutkelp pco2air=400")

Plots.plot!([0:11;].*365 .+ 365/2,-[fluxmean1;fluxmean6],marker = (:circle,2),line=(:solid,:red),label="annual mean integrated flux without kelp")
Plots.plot!([0:10;].*365 .+ 365/2,-[fluxmean2;fluxmean7],marker = (:circle,2),line=(:solid,:black),label="annual mean integrated flux kelp")
Plots.plot!([0:11;].*365 .+ 365/2,-[fluxmean3;fluxmean4],marker = (:circle,2),line=(:solid,:green),label="annual mean integrated flux without kelp DIC=2226")
Plots.plot!([0:11;].*365 .+ 365/2,-[fluxmean5;fluxmean8],marker = (:circle,2),line=(:solid,:blue),label="annual mean integrated flux without kelp pco2air=400")

# -[fluxmean1;fluxmean6]
# -[fluxmean2;fluxmean7]
# -[fluxmean3;fluxmean4]
# -[fluxmean5;fluxmean8]
# (w[2:end].-w[1:end-1])./w[1:end-1]

# Plots.plot(DICforcingsum1-fluxsum1,line=(:solid,:red),ylabel="(mmol/s)",ylims=(-0.5,0.4),label="integrated DIC forcing+flux without kelp")
# Plots.plot!(DICforcingsum2-fluxsum2,line=(:dot,:red),label="integrated DIC forcing+flux kelp")
# Plots.plot!(DICforcingsum3-fluxsum3,line=(:dot,:red),label="integrated DIC forcing+flux kelp DIC=2226")


# Plots.plot(pco2_1,line=(:dot,:red),ylabel="(mmol/s)",legend=:topright,legendfontsize=5,label="integrated flux without kelp") #, right_margin=30Plots.mm
# Plots.plot!(pco2_2,line=(:dot,:black),label="integrated flux kelp")
# Plots.plot!([results3.t/days;results4.t/days],[pco2_3;pco2_4],line=(:dot,:green),label="integrated flux withoutkelp DIC=2226")
# Plots.plot!(pco2_5,line=(:dot,:blue),label="integrated flux withoutkelp pco2air=400")



# savefig("flux.pdf")
# savefig("DICforcing_flux.pdf")
# savefig("DICforcing_flux_withoutkelp.pdf")
# savefig("DICforcing_flux_kelp.pdf")
# savefig("flux_2226.pdf")
# Plots.plot!(DICforcingsum2+fluxsum2,legend=:right)

DIC1=ncread(file1, "DIC")
DIC2=ncread(file2, "DIC")
DIC3=ncread(file3, "DIC")
DIC4=ncread(file4, "DIC")
DIC5=ncread(file5, "DIC")
DIC6=ncread(file6, "DIC")
DIC7=ncread(file7, "DIC")
DIC8=ncread(file8, "DIC")

DIC1sum=sum((DIC1[1,1,:,:].*volume)[depth_index,:],dims=1)[:]
DIC2sum=sum((DIC2[1,1,:,:].*volume)[depth_index,:],dims=1)[:]
DIC3sum=sum((DIC3[1,1,:,:].*volume)[depth_index,:],dims=1)[:]
DIC4sum=sum((DIC4[1,1,:,:].*volume)[depth_index,:],dims=1)[:]
DIC5sum=sum((DIC5[1,1,:,:].*volume)[depth_index,:],dims=1)[:]
DIC6sum=sum((DIC6[1,1,:,:].*volume)[depth_index,:],dims=1)[:]
DIC7sum=sum((DIC7[1,1,:,:].*volume)[depth_index,:],dims=1)[:]
DIC8sum=sum((DIC8[1,1,:,:].*volume)[depth_index,:],dims=1)[:]

DICmean1 = [mean(DIC1sum[(1:365).+x*365]) for x in [0:5;]]
DICmean2 = [mean(DIC2sum[(1:365).+x*365]) for x in [0:5;]]
DICmean3 = [mean(DIC3sum[(1:365).+x*365]) for x in [0:5;]]
DICmean4 = [mean(DIC4sum[(1:365).+x*365]) for x in [0:5;]]
DICmean5 = [mean(DIC5sum[(1:365).+x*365]) for x in [0:5;]]
DICmean6 = [mean(DIC6sum[(1:365).+x*365]) for x in [0:5;]]
DICmean7 = [mean(DIC7sum[(1:365).+x*365]) for x in [0:4;]]
DICmean8 = [mean(DIC8sum[(1:365).+x*365]) for x in [0:5;]]

################################    Integrated DIC gradient 
# itp1=linear_interpolation(1:2191,DIC1sum)
# DICgradient1=only.(Interpolations.gradient.(Ref(itp1),[x for x in 1:2191]))         #https://discourse.julialang.org/t/differentiation-without-explicit-function-np-gradient/57784/2
# Plots.plot(DICgradient1/1day,line=(:dot,:red),ylims=(-0.4,0.3),ylabel="(mmol/s)",label="dDICdt (mmol/s) without kelp, mean="*string(round(mean(DICgradient1/1day),digits=4))*", std="*string(round(std(DICgradient1/1day),digits=4)),legend=:topleft,right_margin=30Plots.mm) # based on tracer conservation equation, has to consider diffusion term 

# itp2=linear_interpolation(1:2191,DIC2sum)
# DICgradient2=only.(Interpolations.gradient.(Ref(itp2),[x for x in 1:2191]))         #https://discourse.julialang.org/t/differentiation-without-explicit-function-np-gradient/57784/2
# Plots.plot!(DICgradient2/1day,line=(:dot,:black),label="dDICdt (mmol/s) kelp, mean="*string(round(mean(DICgradient2/1day),digits=4))*", std="*string(round(std(DICgradient2/1day),digits=4))) 

# itp3=linear_interpolation(1:2191,DIC3sum)
# DICgradient3=only.(Interpolations.gradient.(Ref(itp3),[x for x in 1:2191]))         #https://discourse.julialang.org/t/differentiation-without-explicit-function-np-gradient/57784/2
# Plots.plot!(DICgradient3/1day,line=(:dot,:green),label="dDICdt (mmol/s) without kelp DIC=2226, mean="*string(round(mean(DICgradient3/1day),digits=4))*", std="*string(round(std(DICgradient3/1day),digits=4))) 





# savefig("dDICdt_withoutkelp.pdf")
# savefig("dDICdt_kelp.pdf")
################################

################################    Integrated DIC
# pp=twinx()
# Plots.plot(pp,DIC1sum,ylims=(5.2e8,5.3e8),ylabel="DIC sum (mmol)", line=(:solid,:red),label="integrated DIC without kelp",legendfontsize=5,legend=:bottomright,right_margin=30Plots.mm)
# Plots.plot!(pp,DIC2sum,ylims=(5.2e8,5.3e8),ylabel="DIC sum (mmol)", line=(:solid,:black),label="integrated DIC kelp",legend=:bottomright)
# Plots.plot!(pp,DIC3sum,ylims=(5.20e8,5.38e8),ylabel="DIC sum (mmol)", line=(:solid,:green),label="integrated DIC without kelp DIC=2226",legend=:bottomright)

# Plots.plot([results1.t/days;results6.t/days],[DIC1sum;DIC6sum],ylims=(5.2e8,5.4e8),xlabel="Days",ylabel="DIC sum (mmol)", line=(:solid,:red),label="integrated DIC without kelp",legendfontsize=5,legend=:topright)
# Plots.plot!([results2.t/days;results7.t/days],[DIC2sum;DIC7sum],line=(:solid,:black),label="integrated DIC kelp")
# Plots.plot!([results3.t/days;results4.t/days],[DIC3sum;DIC4sum], label="integrated DIC without kelp DIC=2226",line=(:solid,:green))
# Plots.plot!([results5.t/days;results8.t/days],[DIC5sum;DIC8sum],line=(:solid,:blue),label="integrated DIC without kelp pco2air=400")

# Plots.plot!(moving_avg([DIC1sum;DIC6sum[2:end]]),label="annual mean integrated DIC without kelp")          #moving average
# Plots.plot!(moving_avg([DIC2sum;DIC7sum[2:end]]),label="annual mean integrated DIC kelp")
# Plots.plot!(moving_avg([DIC3sum;DIC4sum[2:end]]),label="annual mean integrated DIC without kelp DIC=2226")
# Plots.plot!(moving_avg([DIC5sum;DIC8sum[2:end]]),label="annual mean integrated DIC without kelp pco2air=400")

# Plots.plot!([0:11;].*365 .+ 365/2,[DICmean1;DICmean6],marker = (:circle,2),line=(:dash,:red),label="annual mean integrated DIC without kelp")
# Plots.plot!([0:10;].*365 .+ 365/2,[DICmean2;DICmean7],marker = (:circle,2),line=(:dash,:black),label="annual mean integrated DIC kelp")
# Plots.plot!([0:11;].*365 .+ 365/2,[DICmean3;DICmean4],marker = (:circle,2),line=(:dash,:green),label="annual mean integrated DIC without kelp DIC=2226")
# Plots.plot!([0:11;].*365 .+ 365/2,[DICmean5;DICmean8],marker = (:circle,2),line=(:dash,:blue),label="annual mean integrated DIC without kelp pco2air=400")


# Plots.plot([results1.t/days;results6.t/days],[DIC1[1,1,Nz,:];DIC6[1,1,Nz,:]],xlabel="Days",ylabel="DIC sum (mmol)", line=(:solid,:red),label="surface DIC without kelp",legendfontsize=5,legend=:topright)
# Plots.plot!([results2.t/days;results7.t/days],[DIC2[1,1,Nz,:];DIC7[1,1,Nz,:]],line=(:solid,:black),label="surface DIC kelp")
# Plots.plot!([results3.t/days;results4.t/days],[DIC3[1,1,Nz,:];DIC4[1,1,Nz,:]], label="surface DIC without kelp DIC=2226",line=(:solid,:green))
# Plots.plot!([results5.t/days;results8.t/days],[DIC5[1,1,Nz,:];DIC8[1,1,Nz,:]],line=(:solid,:blue),label="surface DIC without kelp pco2air=400")



# savefig("dDICdt_withoutkelp.pdf")
# savefig("DIC_dDICdt_both.pdf")
# savefig("DIC_dDICdt_DIC2226.pdf")
# savefig("integrated_DIC.pdf")
# savefig("surface_DIC.pdf")
################################


function NPP(result,output)
    for (iter, t) in enumerate(result.t)
        #tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DIC, :ALK, :DOM, :OXY,:PAR)
        for j in 1:Nz
            NO₃ = result.results[1, 1, 1, j, iter] 
            NH₄ = result.results[2, 1, 1, j, iter] 
            P = result.results[3, 1, 1, j, iter] 
            Z = result.results[4, 1, 1, j, iter] 
            D = result.results[5, 1, 1, j, iter] 
            DD = result.results[6, 1, 1, j, iter] 
            # DIC = result.results[7, 1, 1, j, iter]   
            # ALK = result.results[8, 1, 1, j, iter]
            DOM = result.results[9, 1, 1, j, iter]
            PAR = result.results[11, 1, 1, j, iter]
            # output[j,floor(Int,t/days)+1] = LOBSTER.DIC_forcing(0.5*Lx, 0.5*Ly, result.z[j], t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR, LOBSTER.defaults)
            output[j,iter] = LOBSTER.NPP_forcing(0.5*Lx, 0.5*Ly, result.z[j], t, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults)
        end
            # NPP_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR::AbstractFloat, params) = (1-params.γ)*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P         end
    end
end
NPP_1 = zeros(Nz,floor(Int,duration/days)+1)
NPP(results1,NPP_1)


# pco21 = jldopen(file1pco2, "r")
# pco22 = jldopen(file2pco2, "r") # pco21["pCO₂_record"]
# U_10 = 10.0
# pCO2_air = 413.3
# flux1 = 7.7e-4*U_10^2*(pco21["pCO₂_record"].-pCO2_air)/years*1000   #(365*24*3600/1000) # from mol/m^2/year to mmol/m^2/s    should be multiply by 1000 to get mmol !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# flux2 = 7.7e-4*U_10^2*(pco22["pCO₂_record"].-pCO2_air)/years*1000   #(365*24*3600/1000) # mmol/m^2/s







