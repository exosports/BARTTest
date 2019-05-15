#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import scipy.constants   as const
import scipy.interpolate as si
from opacityconv import bin2np
import voigt

"""
This function generates a Voigt profile based on the Doppler and Lorentzian 
HWHMs used by Transit, and the theoretical Voigt profile for a given molecule 
at some pressure and temperature. These are plotted with the relevant layer in 
the opacity table for comparison. The FWHMs of the profiles are calculated for 
comparison as well.
"""

def comp(molmass=16.04, press=3.3516e-01, temp=1442.58, 
         dia=np.array([3.2e-8, 1.55e-8]), 
         wavenum=4368.350283, 
         atm=np.array([[0.0001, 18.01528], [0.9999, 14.0067]]), 
         opacity='../code-output/01BART/f04broadening/broadening.opt', 
         ipress=35, itemp=9, 
         imol=0, savename='../results/01BART/f04voigt_comp'):
    """
    molmass : float.  Molar mass of the molecule responsible for the spectral 
					  line.
    press   : float.  Pressure of the layer with the molecule.
    temp    : float.  Temperature of the layer with the molecule.
    dia     : float.  Collision diameter of the molecule.
    wavenum : float.  Wavenumber of the line.
    atm     : array.  Contains the abundance and molar mass of each molecule 
                      present in the layer. Shape is (n_molecules, 2)
    opacity : string of path/to/file for the binary opacity file produced by 
                      Transit.
    ipress  : int.    Index of the pressure in the opacity array.
    itemp   : int.    Index of the temperature in the opacity array.
    imol    : int.    Index of the molecule in the opacity array.
	savename: string. Name of produced plot of comparison of profiles.
    """
    # scipy constants
    c  = const.c   * 100 #cm s-1
    k  = const.k   * 1e7
    R  = const.R   * 1e7
    e  = const.e   *  10. * const.c #Coulomb --> statC
    me = const.m_e * 1000.
    h  = const.h   * 1e7
    amu = const.physical_constants['atomic mass constant'][0] * 1000
    Nava = const.N_A

    # Molecule information:
    mass     = molmass * amu        # mass in g
    presscgs = press * 1000000      # pressure in cgs units
    nd       = presscgs/(k*temp)    # number density in cm-3

    # Line information
    wvlen = 10000./wavenum  #wavelength in microns

    # Load Transit opacity file info
    hdr, molID, temps, press, wns, optable = bin2np(opacity, 
                                                    opacity.rsplit('/', 1)[0]+\
                                                    '/opacity')

    # Transit profile
    opa   = optable[ipress,itemp,imol]
    wnhmt = wns[opa > np.amax(opa)/2.0]
    fwhmt = wnhmt[-1] - wnhmt[0]

    # Theoretical Voigt profile calculation (in wavenumber space)
    # Equation for alpha from the Transit Code Manual:
    alpha = wavenum/c * (2.*k*temp*np.log(2)/mass)**0.5
    
    # Equation for gamma from the Transit Code Manual:
    gamma = np.sum((np.pi * np.mean(dia)**2) * (2.*k*temp/np.pi**3/c**2)**0.5 *\
                   (atm[:,0] * nd) *                                        \
                   (1./(molmass*amu) + 1./(atm[:,1]*amu))**0.5)
    #gamma = np.sum((np.pi * (dia)**2) * (2.*k*temp/np.pi**3/c**2)**0.5 *    \
    #               (np.mean(atm[:,0]*atm[:,1]) * atm[:,0] * nd)/atm[:,1] *  \
    #               (1./(molmass*amu) + 1./(atm[:,1]*amu))**0.5)

    # Generate profile
    wnhires = np.linspace(wns[0], wns[-1], 10000)
    prof  = voigt.V(wnhires, alpha, gamma, wavenum)

    Elow = 2183.6851
    gf   = 7.025259715158783e-08
    Z    = 1939.68009857

    prof  = prof * const.pi * e**2 / c**2 / me * nd * atm[0,0] * \
            gf / Z * np.exp(-h*c*Elow/k/temp) * (1 - np.exp(-h*c*wavenum/k/temp))

    # Bin down to `wns`
    resamp = si.interp1d(wnhires, prof, bounds_error=False, fill_value=0)
    prof2  = resamp(wns)

    # Calculate FWHM for funsies
    #wnhm  = wns[prof > np.amax(prof)/2.0]
    #fwhm  = wnhm[-1] - wnhm[0]

    # Make plot
    fig1 = plt.figure(1)
    frame1 = fig1.add_axes((.1, .3, .8, .6))
    plt.plot(wns, prof2, "-", lw=1.5, \
             color="blue", label="Theoretical")
    plt.plot(wns, opa, ":", lw=3.5, \
             color="red", label="Transit")
    #plt.plot(wvrng, profc/np.amax(profc), "--", lw=1.5, color="black", \
    #         label="Calculated")
    plt.xlim(4368.12, 4368.5)
    plt.legend(loc="upper left")
    plt.title("Voigt Profile Comparison")
    plt.xlabel("Wavelength  (um)")
    frame1.set_xticklabels([])
    frame2 = fig1.add_axes((.1, .1, .8, .2))
    #frame2.set_yticklabels([-0.0010, -0.0005, 0.0000, 0.0005, 0.0010])
    plt.plot(wns, prof2 - opa, 'or')
    plt.xlim(4368.12, 4368.5)
    plt.savefig(savename+".png")
    plt.close()


if __name__ == '__main__':
    """
    Run the above code to make the plot
    """
    comp()


