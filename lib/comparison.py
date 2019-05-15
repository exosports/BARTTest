#! /usr/bin/env python

import sys
import numpy as np
import scipy.interpolate as si
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
sys.path.append("../../BART/modules/transit/scripts/")
import readtransit as rt

plt.ion()

def compspec(transit, rhd, geo, atm, outdir=None, titles=False):
    """
    This function produces a plot of Transit and RHD (or other RT code) spectra.

    Inputs
    ------
    transit: string. path/to/file for transit spectrum data.
    rhd    : string. path/to/file for RHD (or other RT code) spectrum data.
    geo    : string. Viewing geometry. 'eclipse' or 'transit'
    atm    : string. Atmosphere's PT profile. 'iso', 'inv', or 'noi'
    outdir : string. path/to/output. Default is execution directory
    titles : bool.   Determines whether to include titles on plots or not.

    Revisions
    ---------
    2017-10-16  mhimes          Initial implementation.
    2018-02-03  raechel         Improved plotting, added residuals subplot.
    2019-04-01  mhimes          Merged into BARTTest.
    """
    # Load transit data
    wlength, flux = rt.readspectrum(transit, 0)

    # Load RHD data
    data  = open(rhd, "r")
    lines = data.readlines()

    if geo == 'eclipse':
        lines = lines[3:]
        
        wlengthb = np.zeros(len(lines), dtype=float)
        fluxb    = np.zeros(len(lines), dtype=float)
        
        for i in range(len(lines)):
            line = lines[i].split()
            wlengthb[i] = float(line[1])
            fluxb[i]    = float(line[2]) * 2.998e10
    elif geo == 'transit':
        lines = lines[3:]
        
        wlengthb = np.zeros(len(lines), dtype=float)
        fluxb    = np.zeros(len(lines), dtype=float)
        
        for i in range(len(lines)):
            line = lines[i].split()
            wlengthb[i] = float(line[1])
            fluxb[i]    = float(line[2]) / 100

    # Resample Transit to RHD's sampling
    rep    = si.splrep(wlength[::-1], flux[::-1]) # Make spline representation
    resamp = si.splev(wlengthb[::-1], rep)        # Evaluate at new x-values
    resamp = resamp[::-1]
    # Set plot title and file output name
    if geo == 'transit':
        if atm == 'inv':
            titlenm = "Inverted Transmission"
            filenm  = "Inverted_transmission_comp.png"
        elif atm == 'iso':
            titlenm = "Isothermal Transmission"
            filenm  = "isothermal_transmission_comp.png"
        elif atm == 'noi':
            titlenm = "Non-inverted Transmission"
            filenm  = "noninverted_transmission_comp.png"
        else:
            print("Wrong `atm` specification. Use 'inv', 'iso', or 'noi'.\n")
    elif geo == 'eclipse':
        if atm == 'inv':
            titlenm = "Inverted Emission"
            filenm  = "Inverted_emission_comp.png"
        elif atm == 'iso':
            titlenm = "Isothermal Emission"
            filenm  = "isothermal_emission_comp.png"
        elif atm == 'noi':
            titlenm = "Non-inverted Emission"
            filenm  = "noninverted_emission_comp.png"
        else:
            print("Wrong `atm` specification. Use 'inv', 'iso', or 'noi'.\n")
    else:
        print("Wrong `geo` specification. Use 'transit' or 'eclipse'.\n")
        sys.exit()
    
    # Plot it
    fig0 = plt.figure(0, (8,5))
    plt.clf()
    frame1 = fig0.add_axes((.14, .3, .8, .65))
    if titles==True:
        plt.title(titlenm)

    if geo == 'transit':
        flux   = 100*flux   # convert to %
        fluxb  = 100*fluxb
        resamp = 100*resamp
    if geo == 'eclipse':
        if atm == 'iso':
            plt.plot(wlength,  flux,   "lightblue", 
                     label='Transit - high resolution', linewidth = 10.0)
            plt.plot(wlengthb, resamp, "royalblue", 
                     label='Transit - binned down to RHD')
            plt.plot(wlengthb, fluxb,  "firebrick", label='RHD')
        else:
            plt.plot(wlength,  flux,   "lightblue", 
                     label='Transit - high resolution')
            plt.plot(wlengthb, resamp, "royalblue", 
                     label='Transit - binned down to RHD')
            plt.plot(wlengthb, fluxb,  "firebrick", label='RHD')
    else:
        plt.plot(wlength,  flux,   "lightblue", 
                 label='Transit - high resolution')
        plt.plot(wlengthb, resamp, "royalblue", 
                 label='Transit - binned down to RHD')
        plt.plot(wlengthb, fluxb,  "firebrick", label='RHD')


    plt.xlabel(u"Wavelength  (\u03bcm)")
    if geo == 'transit':
        plt.ylabel("Modulation  (%)") 
    else:
        plt.ylabel("Flux (erg s$^{-1}$ cm$^{-1}$)")

    #plot black body curves
    h  = 6.62607004e-27     # erg*s
    c  = 2.9979245800e10    # cm/s
    k  = 1.38064852e-16     # erg/K

    # Set min and max temperatures for plotting
    if atm == 'inv':
        T1 =  968.60        # K
        T2 =  1243.05       # K
    elif atm == 'iso':
        T1 =  1100.0               
        T2 =  1100.0
    elif atm == 'noi':
        T1 =  882.93
        T2 =  1408.70

    w   = np.linspace(10000./11.0, 10000./1.000, 10000)   # cm^-1
    
    # Blackbody associated with min/max temperature
    a   = np.pi * 2 * h * (c**2) * (w**3)    
    b1  = (h*c*w) / (k*T1)
    Bb1 = a / (np.exp(b1) - 1.0)
  
    b2  = (h*c*w) / (k*T2)
    Bb2 = a / (np.exp(b2) - 1.0)
    
    l   = 10000. / w
    
    T1s = str(T1)+' K'
    T2s = str(T2)+' K'
    
    if geo == 'eclipse':
        frame1.set_ylim(0, max(Bb2)+5000)
        if T1 != T2:
            plt.plot(l, Bb1, "r", linestyle =':', 
                     label="Blackbody at "+T1s, linewidth = 2.0)
            plt.plot(l, Bb2, "b"        , linestyle =':', 
                     label="Blackbody at "+T2s, linewidth = 2.0)
        else:
            plt.plot(l, Bb1, "k"        , linestyle =':', 
                     label="Blackbody at "+T2s, linewidth = 2.0)
        plt.legend(loc='upper left', prop={'size':8})
    elif geo == 'transit':
        plt.legend(loc='lower right', prop={'size':8})

    frame1.set_xscale('log')
    frame1.set_xlim(1, 11.0)
    frame1.set_xticklabels([])
    frame1.yaxis.set_label_coords(-0.1, 0.5)
    # Set y axis limits
    if geo == 'transit':
        frame1.set_ylim(np.amin(flux)*0.99, np.amax(flux)*1.01)
        frame1.yaxis.set_major_locator(MaxNLocator(nbins='6', prune='lower'))
    else:
        frame1.set_ylim(np.amin(Bb1)*0.9,  np.amax(Bb2)*1.1)


    frame2 = fig0.add_axes((.14, .1, .8, .2))
    
    #residual plots
    resid = np.zeros(len(lines), dtype=float)

    mr = max(resamp)
    mf = max(fluxb)
    m  = max(mr, mf)

    # Residuals, in units of %
    resid = (resamp - fluxb) / np.max(np.concatenate((resamp, fluxb))) * 100

    plt.plot(wlengthb, resid, "k", linestyle = ":")
    plt.ylabel('Residuals (%)')
    plt.xlabel(u"Wavelength  (\u00b5m)")

    frame2.set_xscale('log')
    frame2.set_xlim(0, 11.0)

    frame2.set_ylim(np.amin(resid) - 0.2*np.abs(np.amin(resid)), 
                    np.amax(resid) + 0.2*np.abs(np.amax(resid)))
    
    frame2.set_xticklabels([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    frame2.yaxis.set_major_locator(MaxNLocator(nbins = '5', prune='upper'))
    if(geo == 'eclipse') and (atm == 'iso'):
        frame2.yaxis.get_major_ticks()[-1].label1.set_visible(False)
    frame2.xaxis.set_major_locator(plt.MultipleLocator(1))

    frame2.yaxis.set_label_coords(-0.1, 0.5)

    if(geo == 'eclipse') and (atm == 'iso'):
        frame3 = fig0.add_axes((.57, .38, .25, .25))
        frame3.set_ylim(23000, 24000)
        frame3.set_xscale('log')
        frame3.set_xlim(3.5, 6.5)
        frame3.set_xticklabels([4, 4.5, 5, 5.5, 6])
        frame3.xaxis.set_major_locator(plt.MultipleLocator(1))
        plt.plot(wlength,  flux,   "lightblue", linewidth = 10.0)
        plt.plot(wlengthb, resamp, "royalblue", linewidth = 2.0)
        plt.plot(wlengthb, fluxb,  "firebrick", linewidth = 2.0)
        plt.plot(l, Bb1, "k", linestyle = ':',  linewidth = 4.0)

    if outdir!=None:
        plt.savefig(outdir+filenm)
    else:
        plt.savefig(filenm)
    plt.close()


if __name__ == "__main__":
    """
    Produce the plots using the above function
    """
    compspec('../code-output/01BART/c03hjclearinv/' + \
             'inv_emission_spectrum.dat', 
             '../code-output/02RHD/full_inv_emission.dat', 
             'eclipse', 'inv', 
             '../results/plots/')

    compspec('../code-output/01BART/c01hjcleariso/' + \
             'iso_emission_spectrum.dat', 
             '../code-output/02RHD/full_iso_emission.dat', 
             'eclipse', 'iso', 
             '../results/plots/')

    compspec('../code-output/01BART/c02hjclearnoinv/' + \
             'noinv_emission_spectrum.dat', 
             '../code-output/02RHD/full_noinv_emission.dat', 
             'eclipse', 'noi', 
             '../results/plots/')

    compspec('../code-output/01BART/c03hjclearinv/' + \
             'inv_transmission_spectrum.dat', 
             '../code-output/02RHD/full_inv_transit.dat', 
             'transit', 'inv', 
             '../results/plots/')

    compspec('../code-output/01BART/c01hjcleariso/' + \
             'iso_transmission_spectrum.dat', 
             '../code-output/02RHD/full_iso_transit.dat', 
             'transit', 'iso', 
             '../results/plots/')

    compspec('../code-output/01BART/c02hjclearnoinv/' + \
             'noinv_transmission_spectrum.dat', 
             '../code-output/02RHD/full_noinv_transit.dat', 
             'transit', 'noi', 
             '../results/plots/')

