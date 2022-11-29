
pH=8.0 # initial pH value guess for air-sea flux calculation
U_10 = 10.0  # 10m/s for CO2 flux calculation only
pCO2_air = 413.3 #  in Jan 2020   https://www.co2.earth/  for CO2 flux calculation only
ρₒ = 1026


function air_sea_flux(DIC, ALK, temp,S) # has to be only including x,y,t without z, because this will apply to the z direction. f(x, y, t) on z-boundaries.
    #https://clima.github.io/OceananigansDocumentation/stable/model_setup/boundary_conditions/
    #https://biocycle.atmos.colostate.edu/shiny/carbonate/
    #TEMP =params.TEMP + 273.15 # temperature from Celsius to Kelvin
    TEMP =temp+ 273.15
    #ALK *= 1.e-6 # microequivalents to equivalents
    #DIC *= 1.e-6 # micromoles to moles
    ALK = ALK*1.e-3/ρₒ # microequivalents to equivalents  from mmol/m^-3 to mol/kg
    DIC = DIC*1.e-3/ρₒ # micromoles to moles    

    #S = salinity_itp(mod(t,364days))
    Boron = 1.179e-5*S # Total Boron mole/kg as a fraction of salinity
    
    K0 = exp(-60.2409 + 9345.17/TEMP + 23.3585*log(TEMP/100) + S*(0.023517 - 0.00023656*TEMP + 0.0047036*(TEMP/100)^2))   # mol/kg/atm 
    K1 = exp(2.18867 - 2275.036/TEMP - 1.468591*log(TEMP) + (-0.138681 - 9.33291/TEMP)*sqrt(S) + 0.0726483*S - 0.00574938*S^1.5)
    K2 = exp(-0.84226 - 3741.1288/TEMP -1.437139*log(TEMP) + (-0.128417 - 24.41239/TEMP)*sqrt(S) + 0.1195308*S - 0.0091284*S^1.5)
    KB = exp( (-8966.90 - 2890.51*sqrt(S) - 77.942*S + 1.726*S^1.5 - 0.0993*S^2)/TEMP + (148.0248 + 137.194*sqrt(S) + 1.62247*S) + (-24.4344 - 25.085*sqrt(S) - 0.2474*S)*log(TEMP) + 0.053105*sqrt(S)*TEMP)

    H = 10^(-pH) # initial guess from arg list

    diff_H = H
    tiny_diff_H = 1.e-15
    iter = 0
    CA = 0
    while diff_H > tiny_diff_H && iter < 1000
        H_old = H
        # solve Tans' equation 13 for carbonate alkalinity from TA
        CA = ALK - (KB/(KB + H))*Boron
        # solve quadratic for H (Tans' equation 12)
        a = CA
        b = K1*(CA - DIC)
        c = K1*K2*(CA - 2*DIC)
        H = (-b + sqrt(b^2 - 4*a*c))/(2*a)
        #How different is new estimate from previous one?
        diff_H = abs(H - H_old)
        iter += 1       
    end
    #pH = -log10(H)
    CO2aq = CA/(K1/H + 2*K1*K2/H^2)*1e6 # Eq 11  μmol/kg
    pCO2 = CO2aq/K0 # Eq 4 (converted from atm to ppmv) ppm 
    HCO3 = K0*K1*CO2aq/K0/H       #μmol/kg
    CO3 = K0*K1*K2*CO2aq/K0/H^2   #μmol/kg
    DIC = CO2aq + HCO3 + CO3
    R = DIC/CO3
    #[flux, pCO2, CO2aq, HCO3, CO3, DIC, R]
    flux = 7.7e-4*U_10^2*(pCO2-pCO2_air)/(365*24*3600/1000) # mmol/m^2/s
    return pCO2  #flux      #Wanninkhof 2014 equ.6 positive value means upward flux meaning losing Carbon 
    
end


#air_sea_flux(DIC, ALK, temp,S)
air_sea_flux(2200, 2400, 7,35)
air_sea_flux(2200, 2400, 7,35.5)
air_sea_flux(2200, 2400, 13,35)
air_sea_flux(2200, 2400, 13,35.5)

air_sea_flux(2600, 2400, 7,35)
air_sea_flux(2600, 2400, 7,35.5)
air_sea_flux(2600, 2400, 13,35)
air_sea_flux(2600, 2400, 13,35.5)


air_sea_flux(2200, 2000, 7,35)
air_sea_flux(2200, 2000, 7,35.5)
air_sea_flux(2200, 2000, 13,35)
air_sea_flux(2200, 2000, 13,35.5)

air_sea_flux(2600, 2000, 7,35)
air_sea_flux(2600, 2000, 7,35.5)
air_sea_flux(2600, 2000, 13,35)
air_sea_flux(2600, 2000, 13,35.5)