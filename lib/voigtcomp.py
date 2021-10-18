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

def comp(molmass=18.0105646, press=3.3516e-01, temp=1442.58, 
         dia=np.array([3.2e-8, 1.55e-8]), 
         wavenum=4368.0, 
         atm=np.array([[0.0001, 18.0105646], [0.9999, 14.0067]]),
         ratio=0.9973, 
         Elow=2183.6851, gf=7.026386513565115e-08, Z=2525.51990555, 
         opacity='../code-output/01BART/f04broadening/broadening.opt', 
         ipress=35, itemp=9, imol=0, title=False, 
         savename='../results/01BART/f04voigt_comp', fext='.pdf'):
    """
    molmass : float.  Molar mass of the molecule responsible for the spectral 
					  line.
    press   : float.  Pressure of the layer with the molecule.
    temp    : float.  Temperature of the layer with the molecule.
    dia     : float.  Collision diameter of the molecule.
    wavenum : float.  Wavenumber of the line.
    atm     : array.  Contains the abundance and molar mass of each molecule 
                      present in the layer. Shape is (n_molecules, 2)
    ratio   : float.  Ratio of the isotope.
    Elow    : float.  Energy of the lower state.
    gf      : float.  Weighted oscillator strength.
    Z       : float.  Partition function for the given conditions.
    opacity : string of path/to/file for the binary opacity file produced by 
                      Transit.
    ipress  : int.    Index of the pressure in the opacity array.
    itemp   : int.    Index of the closest temperature in the opacity array.
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

    # Load Transit opacity file info
    hdr, molID, temps, press, wns, optable = bin2np(opacity, 
                                                    opacity.rsplit('/', 1)[0]+\
                                                    '/opacity')

    # Transit profile -- interpolate to the conditions
    opa   = optable[ipress,itemp,imol]
    if round(temp, -2) < temp:
        weight = (temp - round(temp, -2)) / 100
        opa    = (1-weight) * optable[ipress,itemp,imol] + \
                    weight  * optable[ipress,itemp+1,imol]
    else:
        weight = (round(temp, -2) - temp) / 100
        opa    = (1-weight) * optable[ipress,itemp,imol] + \
                    weight  * optable[ipress,itemp-1,imol]

    # Theoretical Voigt profile calculation (in wavenumber space)
    # Equation for Doppler HWHM:
    alpha = wavenum/c * (2.*k*temp*np.log(2)/mass)**0.5
    # Equation for Lorentzian HWHM:
    gamma = np.sum((dia[0]/2. + dia/2.)**2 *                         \
                   (2. * k * temp / np.pi)**0.5 / c *              \
                   nd * atm[:,0] *                              \
                   (1./mass + 1./(atm[:,1] * amu))**0.5)

    # Generate profile
    wnhires = np.linspace(wns[0], wns[-1], 10000)
    prof    = voigt.V(wnhires, alpha, gamma, wavenum)
    # Extinction coefficient
    K = const.pi * e**2 / c**2 / me            *     \
        nd * ratio                             *     \
        gf / Z * np.exp(-h*c*Elow/k/temp)      *     \
        (1 - np.exp(-h*c*wavenum/k/temp))
    # Divide by density -- extinction --> opacity
    prof = K*prof / (presscgs * molmass / k / temp / Nava)
    # Resample to `wns`
    resamp = si.interp1d(wnhires, prof, bounds_error=False, fill_value=0)
    prof2  = resamp(wns)

    # Make plot
    fig1 = plt.figure(1)
    frame1 = fig1.add_axes((.1, .3, .8, .6))
    plt.plot(wns, prof2, "-", lw=1.5, \
             color="blue", label="Theoretical")
    plt.plot(wns, opa, ":", lw=3.5, \
             color="red", label="Transit")
    plt.xlim(4367.5, 4368.5)
    plt.legend(loc="upper left")
    if title:
        plt.title("Voigt Profile Comparison")
    plt.ylabel("Opacity (cm$^2$ g$^{-1}$)")
    frame1.set_xticklabels([])
    yticks = frame1.yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    frame2 = fig1.add_axes((.1, .1, .8, .2))
    plt.plot(wns[opa!=0], ((prof2 - opa)/prof2)[opa!=0], 'or')
    plt.xlim(4367.5, 4368.5)
    plt.xlabel("Wavenumber (cm$^{-1}$)")
    plt.ylabel("Difference (%)")
    yticks = frame2.yaxis.get_major_ticks()
    yticks[-2].label1.set_visible(False)
    plt.savefig(savename+fext, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    """
    Run the above code to make the plot
    """
    comp()


