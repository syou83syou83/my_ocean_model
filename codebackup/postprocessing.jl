
"""
Lobster model based on
2005 A four-dimensional mesoscale map of the spring bloom in the northeast Atlantic (POMME experiment): Results of a prognostic model
2001 Impact of sub-mesoscale physics on production and subduction of phytoplankton in an oligotrophic regime
2012How does dynamical spatial variability impact 234Th-derived estimates of organic export
annual cycle 
"""

using JLD2
using Printf
using Oceananigans
using Plots
using Statistics
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using NetCDF
using HDF5
#using Gnuplot
using Interpolations

#using PyPlot
#using Gadfly
#using PlotlyJS   #interactive

############### Tracer concentration vertical profile

profile = jldopen("profile_subpolar_PARfield_constantlargeKboth.jld2", "r") 
# profile_subpolar_PARfield1e-5, profile_subpolar_PARfield.jld2, profile_subpolar_PARfield_idealizedmld.jld2, profile_subtropical_PARfield_idealizedmld
loc=profile["grid/zᵃᵃᶜ"][4:end-3]
#keys(profile) keys(profile["timeseries"])

iterations = parse.(Int, keys(profile["timeseries/t"])) # turn an array of strings into numbers 
times = [profile["timeseries/t/$iter"] for iter in iterations]
#intro = searchsortedfirst(times, 0)   #23.5hours

#NO₃=profile["timeseries/NO₃"]
#NH₄=profile["timeseries/NH₄"]
#P=profile["timeseries/P"]
#Z=profile["timeseries/Z"]
#D=profile["timeseries/D"]
#DD=profile["timeseries/DD"]
#DOM=profile["timeseries/DOM"]
#DIC=profile["timeseries/DIC"]
#ALK=profile["timeseries/ALK"]

Nz=profile["grid/Nz"]
NO₃_save=zeros(Nz,size(iterations)[1]);
NH₄_save=zeros(Nz,size(iterations)[1]);
P_save=zeros(Nz,size(iterations)[1]);
Z_save=zeros(Nz,size(iterations)[1]);
D_save=zeros(Nz,size(iterations)[1]);
DD_save=zeros(Nz,size(iterations)[1]);
DOM_save=zeros(Nz,size(iterations)[1]);
DIC_save=zeros(Nz,size(iterations)[1]);
ALK_save=zeros(Nz,size(iterations)[1]);
PAR_save=zeros(Nz,size(iterations)[1]);
time_save=zeros(size(iterations)[1])

#anim = @animate for (i, iter) in enumerate(iterations)   #[intro:end])
for (i, iter) in enumerate(iterations)
    @info "Drawing frame $i from iteration $iter..."

    t = profile["timeseries/t/$iter"]
    time_save[i]=t;

    NO₃ = profile["timeseries/NO₃/$iter"][1, 1, :]
    NH₄ = profile["timeseries/NH₄/$iter"][1, 1, :]
    P = profile["timeseries/P/$iter"][1, 1, :]
    Z = profile["timeseries/Z/$iter"][1, 1, :]
    D = profile["timeseries/D/$iter"][1, 1, :]
    DD = profile["timeseries/DD/$iter"][1, 1, :]
    DOM = profile["timeseries/DOM/$iter"][1, 1, :]
    DIC = profile["timeseries/DIC/$iter"][1, 1, :]
    ALK = profile["timeseries/ALK/$iter"][1, 1, :]
    PAR = profile["timeseries/PAR/$iter"][1, 1, :]
    #Budget = NO₃ .+ NH₄ .+ P .+ Z .+ D .+ DD .+ DOM
    #flux = file_flux["timeseries/flux/$iter"][1, 1, Nz]

    NO₃_save[:,i] = NO₃[:]; 
    NH₄_save[:,i] = NH₄[:];
    P_save[:,i] = P[:]; 
    Z_save[:,i] = Z[:]; 
    D_save[:,i] = D[:]; 
    DD_save[:,i] = DD[:]; 
    DOM_save[:,i] = DOM[:]
    DIC_save[:,i] = DIC[:]
    ALK_save[:,i] = ALK[:]
    PAR_save[:,i] = PAR[:]
    #Budget_save[:,i] = Budget[:]
    #flux_save[i] = flux[1]

    NO₃_title = @sprintf("Nitrate , t = %s", prettytime(t)) #"nitrate"
    NH₄_title = "ammonium"
    P_title = "phytoplankton"
    Z_title = "zooplankton"
    D_title = "Detritus"
    DD_title = "Large Detritus"
    DOM_title = "dissolved organic matter"
    DIC_title = "DIC"
    ALK_title = "ALK"
    PAR_title = "PAR"
    #Budget_title = @sprintf("N budget, total = %s", sum(Budget))
    #flux_title = @sprintf("Flux = %s (mol/m^2/y)", flux)
    ## Arrange the plots side-by-side.

"""
    NO₃_plot = Plots.plot(NO₃, loc; ylabel="z(m)")
    NH₄_plot = Plots.plot(NH₄, loc; ylabel="z(m)")
    P_plot = Plots.plot(P, loc; ylabel="z(m)")
    Z_plot = Plots.plot(Z, loc; ylabel="z(m)")
    D_plot = Plots.plot(D, loc; ylabel="z(m)")
    DD_plot = Plots.plot(DD, loc; ylabel="z(m)")
    DOM_plot = Plots.plot(DOM, loc; ylabel="z(m)")
    DIC_plot = Plots.plot(DIC, loc; ylabel="z(m)")
    ALK_plot = Plots.plot(ALK, loc; ylabel="z(m)")
    #flux_plot = plot(time_save[1:i],flux_save[1:i]; ylabel="flux")
    #Budget_plot = plot(Budget, zb; ylabel="z(m)")


    Plots.plot(NO₃_plot, NH₄_plot, P_plot, Z_plot, D_plot, DD_plot, DOM_plot, DIC_plot, ALK_plot,layout=(3, 3), size=(2400, 1000),
         title=[NO₃_title NH₄_title P_title Z_title D_title DD_title DOM_title DIC_title ALK_title]) # S doesn't change much
"""

    #iter == iterations[end] && close(file)
end

#mp4(anim, "annual_subpolar.mp4", fps = 10) # hide


fs=8 #4 front size   #https://stackoverflow.com/questions/57976378/how-to-scale-the-fontsizes-using-plots-jl
kwargs = (xlabel="time (days)", ylabel="z (m)")
NO₃_map=heatmap(time_save/(1day),loc,NO₃_save,titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2),clims=(0,15))
NH₄_map=heatmap(time_save/(1day),loc,NH₄_save,titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2))
P_map=heatmap(time_save/(1day),loc,P_save,titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2)) #,clim=(0.1,1.1)
Z_map=heatmap(time_save/(1day),loc,Z_save,titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2))
D_map=heatmap(time_save/(1day),loc,D_save,titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2))
DD_map=heatmap(time_save/(1day),loc,DD_save,titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2))
DOM_map=heatmap(time_save/(1day),loc,DOM_save,titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2))
DIC_map=heatmap(time_save/(1day),loc,DIC_save,titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2))
ALK_map=heatmap(time_save/(1day),loc,ALK_save,titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2))
PAR_map=heatmap(time_save/(1day),loc,PAR_save,titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2),clims=(0,10))
PAR1_map=heatmap(time_save[1:365]/(1day),loc,log10.(PAR_save[:,1:365]),titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", title="log(PARfield)",xlims=(0,365),clims=(-15,0)) 
savefig("subpolar_logPARfield_constantlargekboth.pdf")

PAR_func = jldopen("subtropical_PARfunc.jld2", "r") # profile_subpolar_PARfield.jld2, profile_subpolar_PARfield_idealizedmld.jld2
heatmap(0:1:364,-PAR_func["depth_chl"][end:-1:1],log10.(PAR_func["PAR"][end:-1:1,:]),title="log(PARfunc)",ylims=(-600,0) ,clims= (-15, 0))   
savefig("subtropical_logPARfunc0.pdf")

#=
using CairoMakie
fig = Figure()
Axis(fig[1, 1], xticks = 0:100:300, yticks = -600:100:0,title = "subpolar_logPARfield")
#ylims!(-500,0)
ax, hm = heatmap(fig[1, 1][1, 1], time_save[1:365]/(1day),loc,transpose(log10.(PAR_save[:,1:365])),colormap = Reverse(:heat),colorrange = (-15, 0))  
Colorbar(fig[1, 1][1, 2], hm)
fig
save("subpolar_logPARfield.pdf", fig) # 

PAR_func = jldopen("subpolar_PARfunc.jld2", "r") # profile_subpolar_PARfield.jld2, profile_subpolar_PARfield_idealizedmld.jld2

fig2 = Figure()
Axis(fig2[1, 1], xticks = 0:100:300, yticks = -600:100:0,title = "subpolar_logPARfunc") #xticks = 0:100:300, yticks = -600:100:0
#ylims!(w,low=-600)
ax, hm = heatmap(fig2[1, 1][1, 1], 0:1:364,-PAR_func["depth_chl"][end:-1:1],transpose(log10.(PAR_func["PAR"][end:-1:1,:])),colormap = Reverse(:heat),colorrange = (-15, 0))   
Colorbar(fig2[1, 1][1, 2], hm)
fig2
=#

save("subtropical_logPARfunc.pdf", fig) # subpolar_logPARfunc.pdf
#jldsave("subpolar_PARfunc.jld2"; PAR, depth_chl)  

#Budget_map=heatmap(time_save/(1day),zb,Budget_save,xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2))
plot(NO₃_map,NH₄_map,P_map,Z_map,D_map,DD_map,DOM_map,DIC_map,ALK_map,title=["Nitrate" "Ammonium" "Phytoplankton" "Zooplankton" "Detritus" "Large Detritus" "DOM" "DIC" "ALK"])
#savefig("annual_cycle_subtropical.pdf")

"""
@gp P_save "w image notit" "set logscale cb" 
@gp :- "set xrange[0:730]" "set yrange[0:150]" 
Gnuplot.save(term="pngcairo size 640,600", output="test.png")
"""
filename_sim = "subpolar_physics.nc"     #subtropical subpolar_physics
mlotst = ncread(filename_sim, "mlotst"); #mixed_layer_depth
mlotst_scale_factor = ncgetatt(filename_sim, "mlotst", "scale_factor") 
mlotst_add_offset = ncgetatt(filename_sim, "mlotst", "add_offset")
mixed_layer_depth = mean(mlotst, dims=(1,2))[1:365]*mlotst_scale_factor.+mlotst_add_offset
#mld_itp = LinearInterpolation(time_series_second, mixed_layer_depth)  
heatmap(time_save/(1day),loc,log10.(P_save),titlefontsize=fs, guidefontsize=fs,tickfontsize=fs,legendfontsize=fs, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2)) #,clim=(0.1,1.1)
plot!(0:729,-[mixed_layer_depth;mixed_layer_depth],label = "MLD",legend=:bottomleft,xlabel="time (days)", ylabel="z (m)",xlims=(0,730),lims=(-600,0))
savefig("subpolar_logP__PARfield_simulated.pdf")

#################### database 


filename_database = "subpolar_BGC.nc"   #subpolar_BGC
Rd_phy = 6.56 # C:N ratio for P molC molN-1
ncinfo(filename_database)
phyc = ncread(filename_database, "phyc");  #phyc
phyc_mean = mean(phyc, dims=(1,2))[1,1,:,:]/Rd_phy #from mmolC/m^3 to mmolN/m^3
depth_phyc = ncread(filename_database, "depth");
heatmap(0:729, -depth_phyc[end:-1:1], phyc_mean[end:-1:1,1:730],interpolate = true)

P_data = heatmap(0:729, -depth_phyc[end:-1:1], log10.(phyc_mean[end:-1:1,1:730]),interpolate = true,xlims=(0,730),ylims=(-600,0),clims=(-2.5,0))
P_data = plot!(0:729,-[mixed_layer_depth;mixed_layer_depth],label = "MLD",legend=:bottomleft,xlabel="time (days)", ylabel="z (m)",xlims=(0,730),ylims=(-600,0))

P_sim=heatmap(time_save/(1day),loc,log10.(P_save),title="Log(P)_sim", xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2),clims=(-2.5,0)) #,clim=(0.1,1.1)
#t_node=[0.,50,86,105,280,364] # subpolar, in days specify piecewiselinear nodes, better include boundary nodes[0.0,364], must be float 
#mld_idealize=[250,420,420,20,20,250] # subpolar, can either extract from Mercator model mld_itp.(t) or specify positive mld at above nodes
t_node=[0.,50,95,100,260,364] # subtropical, 
mld_idealize=[105,135,135,10,10,105] # subtropical,

P_sim=plot!([t_node;t_node.+365],-[mld_idealize;mld_idealize],label = "MLD",legend=:bottomleft,xlabel="time (days)", ylabel="z (m)",xlims=(0,730),ylims=(-600,0))
P_data=plot!([t_node;t_node.+365],-[mld_idealize;mld_idealize],label = "MLD",legend=:bottomleft,xlabel="time (days)", ylabel="z (m)",xlims=(0,730),ylims=(-600,0))

P_sim=plot!(0:729,-[mixed_layer_depth;mixed_layer_depth],label = "MLD",legend=:bottomleft,xlabel="time (days)", ylabel="z (m)",xlims=(0,730),ylims=(-600,0))

plot(P_sim,P_data,layout=(2, 1),size=(1000, 1000),title=["Log(P)_sim" "Log(P)_data"])

#plot(phyc_mean[end:-1:1,1],-depth_phyc[end:-1:1])
savefig("subpolar_logP_PARfield_constantlargeKboth.pdf") 
#annual_cycle_subpolar_logP.pdf, subtropical_logP_PARfield_idealizedmld

savefig("annual_cycle_subtropical_no3_simulated.pdf")

no3 = ncread(filename_database, "no3");  #phyc
no3_mean = mean(no3, dims=(1,2))[1,1,:,:]
#plot(no3_mean[end:-1:1,1],-depth_phyc[end:-1:1])
#no3_mean_itp = interpolate((-depth_phyc[end:-1:1], (1:730)), no3_mean[end:-1:1,1:730],Gridded(Linear()))
#NO₃_data_itp = heatmap(no3_mean_itp[end:-1:1,1:730],ylims=(-600,0),clims=(0,15))

#PAR_itp = interpolate((-depth_chl[end:-1:1], (0:364)day), PAR[end:-1:1,:], Gridded(Linear()))
NO₃_data = heatmap(0:729, -depth_phyc[end:-1:1], no3_mean[end:-1:1,1:730],interpolate = true,titlefontsize=fs,ylims=(-600,0),clims=(0,15))

plot(NO₃_map,NO₃_data,layout=(2, 1),size=(1000, 1000),title=["Nitrate_simulated" "Nitrate_data"])

savefig("annual_cycle_subtropical_no3.pdf")

profile2 = jldopen("pco2_water_subtropical.jld2", "r")
pco2_sim = profile2["pco2_bc"]
spco2 = ncread(filename_database, "spco2");  #phyc
spco2_mean = mean(spco2, dims=(1,2))[1,1,:]

plot(pco2_sim[1,1:731],pco2_sim[2,1:731]/10,label = "pCO2_sim",legend=:bottomright,xlabel="time (days)", ylabel="pCO2 (Pa)",xlims=(0,700),ylims=(0,50))
plot!(spco2_mean[1:731],label = "pCO2_database")
savefig("pco2_sim_subtropical.pdf")
#plot(spco2_mean)
#plot(0:729, -[mixed_layer_depth;mixed_layer_depth],label = "MLD",legend=:bottomleft,xlabel="time (days)", ylabel="z (m)",xlims=(0,730),ylims=(-600,0),right_margin=20Plots.mm)
#savefig("annual_cycle_subpolar_no3_database.pdf")




V_d = -3.47e-5  #  Detritus sedimentation speed   ms⁻¹
V_dd = -50/day  #  Detritus sedimentation speed  -v_dd_min       50m/day=0.0005878  ms⁻¹
molar_mass_C=12 # 12g/mol
G_pump = (D_save[findall(x->x==-98, loc),:]*V_d+DD_save[findall(x->x==-98, loc),:]*V_dd)*Rd_phy*molar_mass_C*1day


###########merge plot of export 
P_map=heatmap(time_save/(1day),loc,log10.(P_save),title="Phytoplankton", legend = :none, xlabel="time (days)", ylabel="z (m)", xlims=(0,365*2),right_margin=20Plots.mm) #,clim=(0.1,1.1)
plot!(0:729, -[mixed_layer_depth;mixed_layer_depth],label = "MLD",legend=:bottomleft,xlabel="time (days)", ylabel="z (m)",xlims=(0,730),ylims=(-600,0),color="blue",right_margin=20Plots.mm)

pp=twinx()
plot!(pp,time_save/(1day), G_pump[:], ylabel="Export [mgC m^-2 d^-1]",label = "Gravitational pump",legend=:bottomright,xlims=(0,730),ylims=(-500,0),color="green")
savefig("Gravitational_pump_and_log10P_MLD_subtropical.pdf")
