#! /usr/bin/env python

import sys
import numpy as np
import scipy.interpolate as si
import matplotlib        as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
sys.path.append("../../BART/modules/transit/scripts/")
import readtransit as rt

plt.ion()

def compspec(specs, codes, geo, atm, outdir=None, titles=False, fext='.pdf'):
    """
    This function produces a plot of supplied spectra.  If multiple are given, 
    it computes the residuals for each spectrum with respect to the average 
    of all supplied spectra.

    Inputs
    ------
    specs  : list, string. path/to/file for spectrum data.
    codes  : list, string. Code names corresponding to `specs`.
    geo    : string. Viewing geometry. 'eclipse' or 'transit'
    atm    : string. Atmosphere's PT profile. 'iso', 'inv', or 'noi'
    outdir : string. path/to/output. Default is execution directory
    titles : bool.   Determines whether to include titles on plots or not.
    fext   : str.    Plot file extension. .pdf or .png

    Revisions
    ---------
    2017-10-16  mhimes          Initial implementation.
    2018-02-03  raechel         Improved plotting, added residuals subplot.
    2019-04-01  mhimes          Merged into BARTTest.
    2021-01-08  mhimes          Refactor to allow more codes to be easily added.
    """
    # Set plot title and file output name
    if geo == 'transit':
        if atm == 'inv':
            titlenm = "Inverted Transmission"
            filenm  = "Inverted_transmission_comp"
        elif atm == 'iso':
            titlenm = "Isothermal Transmission"
            filenm  = "isothermal_transmission_comp"
        elif atm == 'noi':
            titlenm = "Non-inverted Transmission"
            filenm  = "noninverted_transmission_comp"
        else:
            print("Wrong `atm` specification. Use 'inv', 'iso', or 'noi'.\n")
    elif geo == 'eclipse':
        if atm == 'inv':
            titlenm = "Inverted Emission"
            filenm  = "Inverted_emission_comp"
        elif atm == 'iso':
            titlenm = "Isothermal Emission"
            filenm  = "isothermal_emission_comp"
        elif atm == 'noi':
            titlenm = "Non-inverted Emission"
            filenm  = "noninverted_emission_comp"
        else:
            print("Wrong `atm` specification. Use 'inv', 'iso', or 'noi'.\n")
    else:
        raise ValueError("Wrong `geo` specification. Use 'transit' or 'eclipse'.\n")
    
    # Plot it
    fig0 = plt.figure(0, (8,5))
    plt.clf()
    frame1 = fig0.add_axes((.14, .3, .8, .65))
    if titles:
        plt.title(titlenm)

    avgspec = 0 #average spectrum, for residuals
    for i in range(len(specs)):
        wlength, flux = rt.readspectrum(specs[i], 0)
        if geo == 'transit':
            flux *= 100 # convert to %
        avgspec += flux
        plt.plot(wlength, flux, label=codes[i])
    avgspec /= len(specs)

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
        T1 =   911.14       # K
        T2 =  1012.68       # K
    elif atm == 'iso':
        T1 =  1100.0               
        T2 =  1100.0
    elif atm == 'noi':
        T1 =   882.93
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
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    frame1.set_xscale('log')
    frame1.set_xlim(1, 11.0)

    frame1.yaxis.set_label_coords(-0.1, 0.5)

    # For the x-axis to be integers, rather than C * 10^n
    formatter = mpl.ticker.FuncFormatter(lambda y, _: '{:.8g}'.format(y))

    if len(specs) > 1:
        frame1.set_xticklabels([])
        fig0.canvas.draw()
        if frame1.get_yticklabels()[0].get_position()[-1] == 0.0:
            plt.setp(frame1.get_yticklabels()[0], visible=False)
        #residual plots
        frame2 = fig0.add_axes((.14, .1, .8, .2))
        for i in range(len(specs)):
            wlength, flux = rt.readspectrum(specs[i], 0)
            if geo == 'transit':
                flux *= 100 # convert to %
            plt.plot(wlength, flux - avgspec)
        plt.ylabel('Residuals\nw.r.t. average')
        plt.xlabel(u"Wavelength  (\u00b5m)")
        frame2.set_xscale('log')
        frame2.set_xlim(1, 11.0)
        frame2.get_xaxis().set_major_formatter(formatter)
        frame2.get_xaxis().set_minor_formatter(formatter)
        frame2.yaxis.set_label_coords(-0.1, 0.5)
    else:
        frame1.get_xaxis().set_major_formatter(formatter)
        frame1.get_xaxis().set_minor_formatter(formatter)

    if(geo == 'eclipse') and (atm == 'iso'):
        frame3 = fig0.add_axes((.57, .38, .25, .25))
        frame3.set_ylim(23000, 24000)
        frame3.set_xscale('log')
        frame3.set_xlim(3.5, 6.5)
        frame3.get_xaxis().set_major_formatter(formatter)
        frame3.get_xaxis().set_minor_formatter(formatter)
        for i in range(len(specs)):
            wlength, flux = rt.readspectrum(specs[i], 0)
            plt.plot(wlength, flux)
        plt.plot(l, Bb1, "k", linestyle = ':',  linewidth = 4.0)

    if outdir is not None:
        plt.savefig(outdir+filenm+"_"+"_".join(codes)+fext, bbox_inches='tight')
    else:
        plt.savefig(filenm+"_"+"_".join(codes)+fext, bbox_inches='tight')
    plt.close()


def compBarstow(fout, ftransit, fnemesis, fchimera, ftaurex, title=None, fct=1):
    """
    Plots comparisons between Transit, NEMESIS, CHIMERA, and TauREx

    Inputs
    ------
    fout    : string. path/to/output plot file
    ftransit: string. path/to/file for Transit spectrum data.
    fnemesis: string. path/to/file for NEMESIS spectrum data.
    fchimera: string. path/to/file for CHIMERA spectrum data.
    ftaurex : string. path/to/file for TauREx  spectrum data.
    title   : bool.   Determines whether to plot a title
    fct     : int.    Factor to multiply the NEMESIS spectrum by.
    """
    # Load the data
    transit = np.loadtxt(ftransit)
    nemesis = np.loadtxt(fnemesis)
    chimera = np.loadtxt(fchimera)
    taurex  = np.loadtxt(ftaurex)

    # Bin to 0.01um resolution
    if fct==1:
        bins = np.arange(0.5, 10.01, 0.01)
        tr_ibin = np.digitize(transit[:,0], bins)
        tr_bin = np.array([np.mean(100*transit[tr_ibin==i,1]) for i in range(len(bins))])
        ne_ibin = np.digitize(nemesis[:,0], bins)
        ch_ibin = np.digitize(chimera[:,0], bins)
        ta_ibin = np.digitize(taurex [:,0], bins)
        ne_bin = np.array([np.mean(fct*nemesis[ne_ibin==i,1]) for i in range(len(bins))])
        ch_bin = np.array([np.mean(100*chimera[ch_ibin==i,1]) for i in range(len(bins))])
        ta_bin = np.array([np.mean(100*taurex [ta_ibin==i,1]) for i in range(len(bins))])
    else:
        # use nemesis grid for bin midpoints
        midbins = nemesis[:,0][nemesis[:,0] > 0.5]
        lftbins = np.zeros(len(midbins))
        rgtbins = np.zeros(len(midbins))
        lftbins[ 0 ] = midbins[ 0 ] - (midbins[ 1] - midbins[ 0 ])/2.
        lftbins[ 1:] = midbins[ 1:] - (midbins[1:] - midbins[:-1])/2.
        rgtbins[-1 ] = midbins[-1 ] + (midbins[-1] - midbins[-2 ])/2.
        rgtbins[:-1] = midbins[:-1] + (midbins[1:] - midbins[:-1])/2.
        bins    = np.concatenate((lftbins, rgtbins[-1:]))
        tr_ibin = np.digitize(transit[:,0], bins)
        tr_bin  = np.array([np.mean(100*transit[tr_ibin==i,1]) for i in range(len(bins))])
        bins    = bins[:-1]
        tr_bin  = tr_bin[:-1]
        ne_bin  = fct*nemesis[:,1]
        ch_bin  = 100*chimera[:,1]
        ta_bin  = 100*taurex [:,1]

    # Plot
    fig0 = plt.figure(0, (8,5))
    plt.clf()
    frame1 = fig0.add_axes((.14, .3, .8, .65))

    plt.semilogx(bins, tr_bin, label='Transit')
    if fct==1:
        plt.semilogx(bins, ta_bin, label='TauREx')
        plt.semilogx(bins, ne_bin, label='NEMESIS')
        plt.semilogx(bins, ch_bin, label='CHIMERA')
    else:
        plt.semilogx(taurex [:,0], 100*taurex [:,1], label='TauREx')
        plt.semilogx(nemesis[:,0], fct*nemesis[:,1], label='NEMESIS')
        plt.semilogx(chimera[:,0], 100*chimera[:,1], label='CHIMERA')

    if fct==1:
        plt.xlim(0.5, 10.0)
    elif fct==100:
        plt.xlim(1.0, 10.0)
    else:
        raise ValueError("`fct` must be 1 or 100.  Value: "+str(fct))
    frame1.set_xticklabels([])
    plt.ylabel('(Rp/Rs)$^2$ (%)')
    if title is not None:
        plt.title(title)
    plt.legend(loc='best')

    # Residuals
    frame2 = fig0.add_axes((.14, .1, .8, .2))
    if fct==1:
        meanspec = (ta_bin + ne_bin + ch_bin) / 3.
        plt.semilogx(bins, tr_bin - meanspec, label='Transit')
        plt.semilogx(bins, ta_bin - meanspec, label='TauREx')
        plt.semilogx(bins, ne_bin - meanspec, label='NEMESIS')
        plt.semilogx(bins, ch_bin - meanspec, label='CHIMERA')
    else:
        meanspec = (100*taurex[:,1] + fct*nemesis[:,1] + 100*chimera[:,1]) / 3.
        plt.semilogx(bins, tr_bin - meanspec[nemesis[:,0] > 0.5], label='Transit')
        plt.semilogx(taurex [:,0], 100*taurex [:,1] - meanspec, label='TauREx')
        plt.semilogx(nemesis[:,0], fct*nemesis[:,1] - meanspec, label='NEMESIS')
        plt.semilogx(chimera[:,0], 100*chimera[:,1] - meanspec, label='CHIMERA')
    plt.xlim(0.5, 10.0)
    plt.xlabel(u'Wavelength (\u00b5m)')
    plt.ylabel('($\Delta$Rp/Rs)$^2$ (%)')
    formatter = mpl.ticker.FuncFormatter(lambda y, _: '{:.8g}'.format(y))
    frame2.get_xaxis().set_major_formatter(formatter)
    frame2.get_xaxis().set_minor_formatter(formatter)
    plt.savefig(fout, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    """
    Produce the plots using the above functions
    """
    # Add new codes here for comparison plots
    compspec(['../code-output/01BART/c01hjcleariso/iso_uni_emission_spectrum.dat'], 
             ['BART'], 
             'eclipse', 'iso', 
             '../results/plots/')

    compspec(['../code-output/01BART/c02hjclearnoinv/noinv_uni_emission_spectrum.dat'], 
             ['BART'], 
             'eclipse', 'noi', 
             '../results/plots/')

    compspec(['../code-output/01BART/c03hjclearinv/inv_uni_emission_spectrum.dat'], 
             ['BART'], 
             'eclipse', 'inv', 
             '../results/plots/')

    compspec(['../code-output/01BART/c01hjcleariso/iso_transmission_spectrum.dat'], 
             ['BART'], 
             'transit', 'iso', 
             '../results/plots/')

    compspec(['../code-output/01BART/c02hjclearnoinv/noinv_transmission_spectrum.dat'], 
             ['BART'], 
             'transit', 'noi', 
             '../results/plots/')

    compspec(['../code-output/01BART/c03hjclearinv/inv_transmission_spectrum.dat'], 
             ['BART'], 
             'transit', 'inv', 
             '../results/plots/')

    # For Harrington et al. (2021) plots
    compBarstow('../results/01BART/c04_CO_1e-4_1000K.pdf', 
                '../code-output/01BART/c04hjclearisoBarstowEtal/CO_1e-4_1000K_transmission_spectrum.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/single_gas/nemesis_files/nemesis_spectrum_CO_1e-4_1000K_noCIA-noRayleigh.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/single_gas/chimera_files/chimera_spectrum_CO_1e-4_1000K_noCIA-noRayleigh.dat',
                '../../BarstowEtal2020/osfstorage/forward_models/single_gas/taurex_files/taurex_spectrum_CO_1e-4_1000K_noCIA-noRayleigh.dat')
    compBarstow('../results/01BART/c04_CO_1e-4_1500K.pdf', 
                '../code-output/01BART/c04hjclearisoBarstowEtal/CO_1e-4_1500K_transmission_spectrum.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/single_gas/nemesis_files/nemesis_spectrum_CO_1e-4_1500K_noCIA-noRayleigh.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/single_gas/chimera_files/chimera_spectrum_CO_1e-4_1500K_noCIA-noRayleigh.dat',
                '../../BarstowEtal2020/osfstorage/forward_models/single_gas/taurex_files/taurex_spectrum_CO_1e-4_1500K_noCIA-noRayleigh.dat')

    compBarstow('../results/01BART/c04_CO_1e-5_1500K.pdf', 
                '../code-output/01BART/c04hjclearisoBarstowEtal/CO_1e-5_1500K_transmission_spectrum.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/single_gas/nemesis_files/nemesis_spectrum_CO_1e-5_1500K_noCIA-noRayleigh.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/single_gas/chimera_files/chimera_spectrum_CO_1e-5_1500K_noCIA-noRayleigh.dat',
                '../../BarstowEtal2020/osfstorage/forward_models/single_gas/taurex_files/taurex_spectrum_CO_1e-5_1500K_noCIA-noRayleigh.dat')

    compBarstow('../results/01BART/c04_model0.pdf', 
                '../code-output/01BART/c04hjclearisoBarstowEtal/model0_transmission_spectrum.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/realistic/nemesis_hd189733b-cloudfree_jwst_60.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/realistic/chimera_hd189733b-cloudfree_jwst_60.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/realistic/taurex_hd189733b-cloudfree_jwst_60.dat', 
                fct=100)

    compBarstow('../results/01BART/c05_model1.pdf', 
                '../code-output/01BART/c05hjcloudisoBarstowEtal/model1_transmission_spectrum.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/realistic/nemesis_hd189733b-clouds_jwst_60.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/realistic/chimera_hd189733b-clouds_jwst_60.dat', 
                '../../BarstowEtal2020/osfstorage/forward_models/realistic/taurex_hd189733b-clouds_jwst_60.dat', 
                fct=100)






