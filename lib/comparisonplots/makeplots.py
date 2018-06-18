import RHDcomp

# Eclipse comparisons
RHDcomp.compspec('/home/mhimes/retrievaltests/PT-inv-run/PTinv_07oct2017_spectrum.dat', '/home/mhimes/Downloads/RHD_results/inv_ds_emission_noTiso.dat', 'eclipse', 'inv', '/home/mhimes/BARTTest_clean/RHDcomp/')

RHDcomp.compspec('/home/mhimes/retrievaltests/PT-iso-run/PTiso_07oct2017_spectrum.dat', '/home/mhimes/Downloads/RHD_results/iso_ds_emission_noTiso.dat', 'eclipse', 'iso', '/home/mhimes/BARTTest_clean/RHDcomp/')

RHDcomp.compspec('/home/mhimes/retrievaltests/PT-noinv-run/PTnoinv_07oct2017_spectrum.dat', '/home/mhimes/Downloads/RHD_results/noinv_ds_emission_noTiso.dat', 'eclipse', 'noi', '/home/mhimes/BARTTest_clean/RHDcomp/')

# Transit comparisons
RHDcomp.compspec('/home/mhimes/retrievaltests/PT-inv-run-t/PTinv_07oct2017_spectrum.dat', '/home/mhimes/Downloads/RHD_results/full_inv_transit.dat', 'transit', 'inv', '/home/mhimes/BARTTest_clean/RHDcomp/')

RHDcomp.compspec('/home/mhimes/retrievaltests/PT-iso-run-t/PTiso_07oct2017_spectrum.dat', '/home/mhimes/Downloads/RHD_results/full_iso_transit.dat', 'transit', 'iso', '/home/mhimes/BARTTest_clean/RHDcomp/')

RHDcomp.compspec('/home/mhimes/retrievaltests/PT-noinv-run-t/PTnoinv_07oct2017_spectrum.dat', '/home/mhimes/Downloads/RHD_results/full_noinv_transit.dat', 'transit', 'noi', '/home/mhimes/BARTTest_clean/RHDcomp/')

