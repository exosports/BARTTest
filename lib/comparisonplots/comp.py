import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("../../../BART/modules/transit/scripts/")
import readtransit as rt

plt.ion()

def compspec(code1, code2, geo, atm, outdir):
    """
    This function produces a plot comparing the spectra produced by two codes.
    """
    # Load code1 data
    wlength1, flux1 = rt.readspectrum(code1, 0)

    # Load code2 data
    wlength2, flux2 = rt.readspectrum(code2, 0)

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
    
    # Plot it
    plt.figure(0, (8,5))
    plt.clf()
    plt.title(titlenm)
    plt.plot(wlength1, flux1, "b", label='code1')
    plt.plot(wlength2, flux2, "g", label='code2')
    plt.xlabel("Wavelength  (um)")
    if geo == 'transit':
        plt.ylabel("Modulation  (Rp/Rs)^2")
    elif geo == 'eclipse':
        plt.ylabel("Flux (erg s$^-1$ cm$^-1$)")
    else:
        print("No y-label on the plot. Look at the other error message above.\n")
    plt.legend(loc='upper right')
    plt.savefig(outdir+filenm)
    plt.close()

