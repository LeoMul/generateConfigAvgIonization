import numpy as np 

def convolMaxwell(
                  energyGrid, 
                  omega,
                  temperatureGrid):
    
    '''
    This function calculates effective 'function'. Mostly collision strengths.
    
    energyGrid - array of energies (eV) (electron energies post excitation/ionization)
    omega - function to be convolved.
    right now - assume that we're relative to upper level. 
    temperatureGrid - temperatures in eV 
    
    returns ups = int_0 ^infty omega  e^{E/kT} dE /kT
    
    '''
    
    oneOverT = np.power(temperatureGrid,-1) 
    #print(oneOverT)
    ups = np.zeros_like(temperatureGrid)
    
    #ot is in eV-1 - i,e = 1/ kT(eV) 
    for (i,ot) in enumerate(oneOverT):
        #energy is the electron after ionization, so
        #no need to play with the excitation energy. 
        
        maxwellGrid = np.exp( - energyGrid * ot)
        integrand = maxwellGrid * omega 
        
        #call an optimized library
        ups[i] = np.trapz(integrand, energyGrid) * ot    
        
        #this is reserved for the rate calculation.
        #ups[i] *= np.exp(-ionization * ot )

    return ups


def processCSA(final_electron_energy_ev,
               csa_mb,
               temperature_grid_kelvin,
               transition_energy_eV):
    
    #Takes the final (post collision) electron energy(eV)
    #and a cross section in meba barns 
    #and outputs a rate in cm3 s-1
    
    EV_KELVIN = 11604.52
    BOLTZ_SI = 1.380649E-23
    PI = 3.14159265359
    MASS_ELECTRON  = 9.109383E-31
    #a_0^2 = 27.99e-22 m^2 = 27.99 Mb
    
    #possibly missing a factor of pi
    #omega = csa_au * energy_ryd * weight 
    
    #this takes care of any issues with types in the json.
    temperature_grid_ev = np.ones(len(temperature_grid_kelvin))

    for (i,temp) in enumerate(temperature_grid_kelvin):        
        temperature_grid_ev[i] = float(temp) / EV_KELVIN
    
    #In the integration of a CROSS SECTION over energy, there is an implicit 
    #conversion to collision strength via the product \sigma * \epsilon.
    #However, this \epsilon is the INITIAL free electron energy.
    #Therefore we must add the transition energy in before the integation.
    
    #convolve, it will have units of eV*Mb
    initial_energy = final_electron_energy_ev+transition_energy_eV
    
    convol_eV_Mb = convolMaxwell(
                    final_electron_energy_ev,
                    csa_mb*initial_energy,
                    temperature_grid_ev
                                )
    
    #remove ev, convert to cm.
    convol_cm2 = convol_eV_Mb*1e-18 / temperature_grid_ev
    
    #constants and so forth.    
    ff = BOLTZ_SI * 8.0 / (PI * MASS_ELECTRON)
    #m per s
    rate_factor = np.sqrt ( ff * temperature_grid_ev * EV_KELVIN )
    #print(rate_factor)
    #cm per sec
    rate_factor*= 100
    
    #in the integration we did relative to the upper level - so we must incude the contribution
    #from the exponential down here.
    rate = rate_factor * convol_cm2 * np.exp ( -transition_energy_eV / temperature_grid_ev)
    
    
    return rate