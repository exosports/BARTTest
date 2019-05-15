.PHONY: all some tests_noiso test_plusiso oneline fewline multiline abundance isothermal broadening blending multicia plots fin clean

all: bart hitran_linelists oneline fewline multiline broadening abundance blending multicia isothermal enregycons comparisontests plots fin

some: bart oneline fewline multiline broadening abundance blending multicia plots fin

quicktests: oneline fewline multiline broadening abundance blending multicia plots fin

forwardtests: oneline fewline multiline broadening abundance blending multicia isothermal energycons plots fin

comparisontests: comparison_tli comparison_iso comparison_noinv comparison_inv comparison_plots fin

synth_retrievals: retrieval_tli retrieval_iso_e retrieval_iso_t retrieval_noinv_e retrieval_noinv_t retrieval_inv_e retrieval_inv_t fin

hd189: hd189_tli hd189_retrieval fin

bart:
	@echo "\nCloning BART..."
	@if [ ! -d "../BART" ]; then                                              \
		git clone --recursive https://github.com/exosports/BART ../BART/;     \
		echo "Finished cloning BART to a directory parallel to BARTTest.\n";  \
		echo "Compiling BART...";                                             \
		cd ../BART/modules/transit/ && make;                                  \
		cd ../BART/modules/MCcubed/ && make;                                  \
		echo "Finished compiling BART.\n";                                    \
	else                                                                      \
		echo "BART already exists.\n";                                        \
	fi

hitran_linelists:
	@echo "Downloading HITRAN line lists...\n"
	@cd tests/00inputs/par/                                                 &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITEMP-CO2.txt    &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITEMP-CO.txt     &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITEMP-H2O.txt    &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITRAN-CH4.txt    &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITRAN-NH3.txt    &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITRAN-H2.txt
	@echo "Extracting archives...\n"
	@cd tests/00inputs/par/                                                 &&\
	unzip '01_*HITEMP2010.zip'                                              &&\
	unzip '02_*HITEMP2010.zip'                                              &&\
	unzip '05_*HITEMP2010.zip'                                              &&\
	unzip '06_hit12.zip'                                                    &&\
	unzip '11_hit12.zip'                                                    &&\
	unzip '45_hit12.zip'                                                    &&\
	rm -f *.zip
	@echo "Finished retrieving HITRAN line lists.\n"

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
	                                             energycons_3_emission.trc  &&\
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
	../../../BART/BART.py -c iso_emission.brt --justPlots

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


hd189_tli:
	@echo "Running retrieval, HD 189733b: \n"
	@cd tests/r01hd189733b/                                                 &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	pyline.plc

hd189_retrieval:
	@echo "Running retrieval, HD 189733b: \n"
	@cd tests/r01hd189733b/                                                 &&\
	../../../BART/BART.py -c HD189733b.brt -- justOpacity

plots:
	@echo "Making plots..."
	@cd lib/ && makeplots.py
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


