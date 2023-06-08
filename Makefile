.PHONY: all quicktests forwardtests comparisontests synth_retrievals hd189 oneline fewline multiline abundance broadening blending multicia isothermal energycons plots fin clean

all: bart linelists oneline fewline multiline broadening abundance blending multicia isothermal energycons comparisontests plots fin

quicktests: oneline fewline multiline broadening abundance blending multicia plots fin

forwardtests: oneline fewline multiline broadening abundance blending multicia isothermal energycons plots fin

comparisontests: comparison_tli comparison_iso comparison_noinv comparison_inv comparison_BarstowEtal BarstowEtal_OSF comparison_plots fin

synthretrievals: retrieval_tli retrieval_iso_e retrieval_iso_t retrieval_noinv_e retrieval_noinv_t retrieval_inv_e retrieval_inv_t retrieval_BarstowEtal_NEMESIS_clear retrieval_BarstowEtal_NEMESIS_cloud retrieval_BarstowEtal_CHIMERA_clear retrieval_BarstowEtal_CHIMERA_cloud retrieval_BarstowEtal_TauREx_clear retrieval_BarstowEtal_TauREx_cloud fin

BarstowEtalretrievals: retrieval_tli retrieval_BarstowEtal_NEMESIS_clear retrieval_BarstowEtal_NEMESIS_cloud retrieval_BarstowEtal_CHIMERA_clear retrieval_BarstowEtal_CHIMERA_cloud retrieval_BarstowEtal_TauREx_clear retrieval_BarstowEtal_TauREx_cloud fin

hd189: hd189_tli hd189_retrieval fin

bart:
	@if [ ! -d "../BART" ]; then                                              \
		@echo "\nCloning BART...";                                            \
		git clone --recursive https://github.com/exosports/BART ../BART/;     \
		@cd ../BART && git checkout a08ee09c67fc7c571efe51b55caba104720e0db2; \
		echo "Finished cloning BART to a directory parallel to BARTTest.\n";  \
	else                                                                      \
		echo "BART already exists in a directory parallel to BARTTest.\n";    \
	fi
	@echo "Compiling BART..."
	@cd ../BART/modules/transit && make
	@cd ../BART/modules/MCcubed && make
	@echo "Finished compiling BART.\n"


linelists:
	@echo "Downloading line lists...\n"
	@cd tests/00inputs/par/                                                 &&\
	wget https://hitran.org/hitemp/data/bzip2format/05_HITEMP2019.par.bz2   &&\
	wget https://hitran.org/hitemp/data/bzip2format/06_HITEMP2020.par.bz2   &&\
	wget https://zenodo.org/record/3768504/files/CO2_exomol_ucl4000_0.5-500.0um_100-3500K_threshold_0.01_lbl.dat             &&\
	wget https://zenodo.org/record/3768504/files/H2O_exomol_pokazatel_0.24-500.0um_100-3500K_threshold_0.01_lbl.dat
	@echo "Extracting archives...\n"
	@cd tests/00inputs/par/                                                 &&\
	wget -i wget-list_HITEMP-H2O.txt                                        &&\
	wget -i wget-list_HITEMP-CO2.txt                                        &&\
	wget -i wget-list_HITEMP-CO.txt                                         &&\
	wget -i wget-list_HITRAN-CH4.txt                                        &&\
	wget -i wget-list_HITRAN-NH3.txt                                        &&\
	wget -i wget-list_HITRAN-H2.txt                                         &&\
	unzip '01_*HITEMP2010.zip'                                              &&\
	unzip '02_*HITEMP2010.zip'                                              &&\
	unzip '05_*HITEMP2010.zip'                                              &&\
	unzip '06_hit12.zip'                                                    &&\
	unzip '11_hit12.zip'                                                    &&\
	unzip '45_hit12.zip'                                                    &&\
	bunzip2 05_HITEMP2019.par.bz2                                           &&\
	bunzip2 06_HITEMP2020.par.bz2                                           &&\
	rm -f *.zip
	@echo "Finished downloading line lists.\n"

oneline:
	@echo "Running oneline test...\n"
	@cd tests/f01oneline/                                                   &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   oneline.plc                                                          &&\
	../../../BART/modules/transit/transit/transit -c oneline_emission.trc
	@echo "oneline test complete.\n"

fewline:
	@echo "Running fewline test...\n"
	@cd tests/f02fewline/                                                   &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   fewline.plc                                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   ./fewline_emission.trc                                               &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   ./fewline_transmission.trc 
	@echo "fewline test complete.\n"

multiline:
	@echo "Running multiline test...\n"
	@cd tests/f03multiline/                                                 &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   multiline.plc                                                        &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   multiline_emission.trc                                               &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   multiline_transmission.trc
	@echo "multiline test complete.\n"

broadening:
	@echo "Running broadening test...\n"
	@cd tests/f04broadening/                                                &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   broadening.plc                                                       &&\
	../../../BART/modules/transit/transit/transit -c broadening_emission.trc
	@cd lib/ && ./voigtcomp.py
	@echo "broadening test complete.\n"

abundance:
	@echo "Running abundance test...\n"
	@cd tests/f05abundance/                                                 &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   abundance.plc                                                        &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_0_emission.trc                                             &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_1e-4_emission.trc                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_2e-4_emission.trc                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_3e-4_emission.trc                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_4e-4_emission.trc                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_5e-4_emission.trc                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_6e-4_emission.trc                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_7e-4_emission.trc                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_8e-4_emission.trc                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_9e-4_emission.trc                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_1e-3_emission.trc
	@cd ./lib/ && ./abuncomp.py
	@echo "abundance test complete.\n"

blending:
	@echo "Running blending test...\n"
	@cd tests/f06blending/                                                  &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   blending.plc                                                         &&\
	../../../BART/modules/transit/transit/transit -c blending_emission.trc
	@echo "blending test complete.\n"

multicia:
	@echo "Running multicia test...\n"
	@cd tests/f07multicia/                                                  &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   multicia.plc                                                         &&\
	../../../BART/modules/transit/transit/transit -c noCIA_emission.trc     &&\
	../../../BART/modules/transit/transit/transit -c oneCIA_emission.trc    &&\
	../../../BART/modules/transit/transit/transit -c twoCIA_emission.trc
	@echo "multicia test complete.\n"

isothermal:
	@echo "Running isothermal test...\n"
	@if [ ! -f ./tests/00inputs/TLI/CH4_CO_CO2_H2O_full_1-8um.tli ]; then     \
	    cd tests/f08isothermal/                                             &&\
	    ../../../BART/modules/transit/pylineread/src/pylineread.py -c         \
	                                                 isothermal.plc;          \
	fi
	@cd tests/f08isothermal/                                                &&\
	../../../BART/modules/transit/transit/transit -c isothermal_emission.trc
	@echo "isothermal test complete. \n"

energycons:
	@echo "Running energy conservation test...\n"
	@if [ ! -f ./tests/00inputs/TLI/CH4_CO_CO2_H2O_full_1-8um.tli ]; then     \
	    cd tests/f09energycons/                                             &&\
	    ../../../BART/modules/transit/pylineread/src/pylineread.py -c         \
	                                                 energycons.plc;          \
	fi
	@cd tests/f09energycons/                                                &&\
	../../../BART/modules/transit/transit/transit -c                          \
	                                             energycons_1_emission.trc  &&\
	../../../BART/modules/transit/transit/transit -c                          \
	                                             energycons_2_emission.trc  &&\
	../../../BART/modules/transit/transit/transit -c                          \
	                                             energycons_5_emission.trc  &&\
	../../../BART/modules/transit/transit/transit -c                          \
	                                             energycons_10_emission.trc
	@cd ./lib/ && ./energycons.py
	@echo "energy conservation test complete. \n"

comparison_tli:
	@echo "Generating TLI file for comparison tests...\n"
	@if [ ! -f ./tests/00inputs/TLI/CH4_CO_CO2_H2O_NH3_H2_1-11um.tli ]; then  \
	    cd tests/c01hjcleariso/                                             &&\
	    ../../../BART/modules/transit/pylineread/src/pylineread.py -c         \
	                                                 comparison.plc;          \
		echo "TLI file generated for comparison tests.\n";                    \
	else                                                                      \
		echo "TLI file already exists.\n";                                    \
	fi

comparison_iso:
	@echo "Running comparison test, isothermal atmosphere: \n"
	@cd tests/c01hjcleariso/                                                &&\
	../../../BART/modules/transit/transit/transit -c iso_emission.trc       &&\
	../../../BART/modules/transit/transit/transit -c iso_transmission.trc
	@echo "Isothermal comparison test complete. \n"

comparison_noinv:
	@echo "Running comparison test, noninverted atmosphere: \n"
	@cd tests/c02hjclearnoinv/                                                &&\
	../../../BART/modules/transit/transit/transit -c noinv_emission.trc     &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   noinv_transmission.trc
	@echo "Noninverted comparison test complete. \n"

comparison_inv:
	@echo "Running comparison test, inverted atmosphere: \n"
	@cd tests/c03hjclearinv/                                                &&\
	../../../BART/modules/transit/transit/transit -c inv_emission.trc       &&\
	../../../BART/modules/transit/transit/transit -c inv_transmission.trc
	@echo "Inverted comparison test complete. \n"

comparison_BarstowEtal:
	@echo "Running comparison test, Barstow et al. (2020) cases: \n"
	@cd tests/c04hjclearisoBarstowEtal/                                     &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   pyline_BarstowEtal_CO.plc                                            &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   pyline_BarstowEtal_model0_H2O_CO.plc                                 &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   BarstowEtal_CO_1e-4_1000K.trc                                        &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   BarstowEtal_CO_1e-4_1500K.trc                                        &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   BarstowEtal_CO_1e-5_1500K.trc                                        &&\
	../../../BART/modules/transit/transit/transit -c BarstowEtal_model0.trc
	@cd tests/c05hjcloudisoBarstowEtal/                                     &&\
	../../../BART/modules/transit/transit/transit -c BarstowEtal_model1.trc
	@echo "Barstow et al. (2020) comparison test cases complete. \n"

BarstowEtal_OSF:
	@if [ ! -d "../BarstowEtal2020" ]; then                                   \
		@echo "\nCloning Barstow et al. (2020) OSF repo..."                   \
		osf -p 5hg6y clone .;                                                 \
		echo "Finished cloning Barstow et al. (2020) OSF repo to a directory parallel to BARTTest.\n";  \
	else                                                                      \
		echo "Barstow et al. (2020) OSF repo already exists in a directory parallel to BARTTest.\n";    \
	fi

comparison_plots:
	@echo "Making plots for comparison test...\n"
	@cd ./lib/ && ./comparison.py

retrieval_tli:
	@echo "Generating TLI file for synthetic retrieval tests...\n"
	@if [ ! -f ./tests/00inputs/TLI/CH4_CO_CO2_H2O_NH3_H2_1-11um.tli ]; then  \
	    cd tests/s01hjcleariso/                                             &&\
	    ../../../BART/modules/transit/pylineread/src/pylineread.py -c         \
	                                                 retrievals.plc;          \
		echo "TLI file generated for synthetic retrieval tests.\n";           \
	else                                                                      \
		echo "TLI file already exists.\n";                                    \
	fi

retrieval_iso_e:
	@echo "Running retrieval, isothermal atmosphere, eclipse: \n"
	@cd tests/s01hjcleariso/                                                &&\
	../../../BART/BART.py -c iso_emission.brt

retrieval_iso_t:
	@echo "Running retrieval, isothermal atmosphere, transit: \n"
	@cd tests/s01hjcleariso/                                                &&\
	../../../BART/BART.py -c iso_transmission.brt

retrieval_noinv_e:
	@echo "Running retrieval, noninverted atmosphere: \n"
	@cd tests/s02hjclearnoinv/                                              &&\
	../../../BART/BART.py -c noinv_emission.brt

retrieval_noinv_t:
	@echo "Running retrieval, noninverted atmosphere: \n"
	@cd tests/s02hjclearnoinv/                                              &&\
	../../../BART/BART.py -c noinv_transmission.brt

retrieval_inv_e:
	@echo "Running retrieval, inverted atmosphere: \n"
	@cd tests/s03hjclearinv/                                                &&\
	../../../BART/BART.py -c inv_emission.brt

retrieval_inv_t:
	@echo "Running retrieval, inverted atmosphere: \n"
	@cd tests/s03hjclearinv/                                                &&\
	../../../BART/BART.py -c inv_transmission.brt

retrieval_BarstowEtal_NEMESIS_clear:
	@echo "Running retrieval, Barstow et al. (2020) NEMESIS Model 0: \n"
	@cd tests/s04hjclearisoBarstowEtal/                                     &&\
	../../../BART/BART.py -c BarstowEtal_model0_nemesis.brt

retrieval_BarstowEtal_NEMESIS_cloud:
	@echo "Running retrieval, Barstow et al. (2020) NEMESIS Model 1: \n"
	@cd tests/s05hjcloudisoBarstowEtal/                                     &&\
	../../../BART/BART.py -c BarstowEtal_model1_nemesis.brt

retrieval_BarstowEtal_CHIMERA_clear:
	@echo "Running retrieval, Barstow et al. (2020) CHIMERA Model 0: \n"
	@cd tests/s04hjclearisoBarstowEtal/                                     &&\
	../../../BART/BART.py -c BarstowEtal_model0_chimera.brt

retrieval_BarstowEtal_CHIMERA_cloud:
	@echo "Running retrieval, Barstow et al. (2020) CHIMERA Model 1: \n"
	@cd tests/s05hjcloudisoBarstowEtal/                                     &&\
	../../../BART/BART.py -c BarstowEtal_model1_chimera.brt

retrieval_BarstowEtal_TauREx_clear:
	@echo "Running retrieval, Barstow et al. (2020) Tau-REx Model 0: \n"
	@cd tests/s04hjclearisoBarstowEtal/                                     &&\
	../../../BART/BART.py -c BarstowEtal_model0_taurex.brt

retrieval_BarstowEtal_TauREx_cloud:
	@echo "Running retrieval, Barstow et al. (2020) Tau-REx Model 1: \n"
	@cd tests/s05hjcloudisoBarstowEtal/                                     &&\
	../../../BART/BART.py -c BarstowEtal_model1_taurex.brt

hd189_tli:
	@echo "Running retrieval, HD 189733 b: \n"
	@cd tests/r01hd189733b/                                                 &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	HD189733b.plc

hd189_opacity:
	@echo "Generating opacity table for HD 189733 b retrieval: \n"
	@cd tests/r01hd189733b/                                                 &&\
	../../../BART/BART.py -c HD189733b.brt --justOpacity

hd189_retrieval:
	@echo "Running retrieval, HD 189733 b: \n"
	@cd tests/r01hd189733b/                                                 &&\
	../../../BART/BART.py -c HD189733b.brt

plots:
	@echo "Making plots..."
	@cd lib/ && ./makeplots.py
	@echo "Plotting complete.\n"

retrievalplots:
	@echo "Making synthetic retrieval plots..."
	@cd lib/ && ./retrievalplots.py && ./retrievalplots_BarstowEtal.py
	@echo "Plotting complete.\n"

fin:
	@echo "Done!\n"

clean:
	@echo "Deleting files created from running BARTTest.\n"
	@echo "Deleting .pyc files..."
	@cd lib/                            &&\
	rm -f *.pyc
	@echo "Deleting BART output files..."
	@cd code-output/01BART/             &&\
	rm -rf ./*
	@echo "Deleting BART results..."
	@cd results/01BART/                 &&\
	rm -rf ./*
	@echo "Deleting downloaded line list files..."
	@cd tests/00inputs/par/             &&\
	rm -f 01* 02* 05* 06* 11* 45*
	@echo "Deleting TLI files...\n"
	@cd tests/00inputs/TLI/             &&\
	rm -f *.tli
	@echo "BARTTest is now back to its base state."


