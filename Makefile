.PHONY: all some tests_noiso test_plusiso oneline fewline multiline abundance isothermal broadening blending multicia plots fin clean

all: bart hitran_linelists oneline fewline multiline broadening abundance blending multicia isothermal comparison plots fin

some: bart oneline fewline multiline broadening abundance blending multicia plots fin

tests_noiso: oneline fewline multiline broadening abundance blending multicia plots fin

tests_plusiso: oneline fewline multiline broadening abundance blending multicia isothermal plots fin

bart:
	@echo "\nCloning BART..."
	git clone --recursive https://github.com/exosports/BART ../BART/
	@echo "Finished cloning BART to a directory parallel to BARTTest.\n"
	@echo "Compiling BART..."
	@cd ../BART/modules/transit/ && make
	@cd ../BART/modules/MCcubed/ && make
	@echo "Finished compiling BART.\n"

hitran_linelists:
	@echo "Downloading HITRAN line lists...\n"
	@cd tests/00inputs/par/                                                 &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITEMP-CO2.txt    &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITEMP-CO.txt     &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITEMP-H2O.txt    &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITRAN-CH4.txt
	@echo "Extracting archives...\n"
	@cd tests/00inputs/par/                                                 &&\
	unzip '01_*HITEMP2010.zip'                                              &&\
	unzip '02_*HITEMP2010.zip'                                              &&\
	unzip '05_*HITEMP2010.zip'                                              &&\
	unzip '06_hit08.zip'                                                    &&\
	@rm -f *.zip
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
	@cd lib/ && python -c 'import voigtcomp;voigtcomp.comp()'
	@echo "broadening test complete.\n"

abundance:
	@echo "Running abundance test...\n"
	@cd tests/f05abundance/                                                 &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   abundance_line.plc                                                   &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   abundance_noline.plc                                                 &&\
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
	   abundance_1e-3_emission.trc                                          &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   abundance_noline_emission.trc
	@cd lib/ && comparison.py
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
	@cd tests/f08isothermal/                                                &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   isothermal.plc                                                       &&\
	../../../BART/modules/transit/transit/transit -c isothermal_emission.trc
	@echo "isothermal test complete. \n"

comparison:
	@echo "Running comparison test cases...\n"
	@cd tests/f09comparison/                                                &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   comparison.plc                                                       &&\
	../../../BART/modules/transit/transit/transit -c iso_emission.trc       &&\
	../../../BART/modules/transit/transit/transit -c iso_transmission.trc   &&\
	../../../BART/modules/transit/transit/transit -c noinv_emission.trc     &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   noinv_transmission.trc                                               &&\
	../../../BART/modules/transit/transit/transit -c inv_emission.trc       &&\
	../../../BART/modules/transit/transit/transit -c inv_transmission.trc
	@echo "Comparison test cases complete. \n"

comparison_tli:
	@echo "Generating TLI file for comparison tests...\n"
	@cd tests/f09comparison                                                 &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   comparison.plc
	@echo "TLI file generated for comparison tests.\n"

comparison_iso:
	@echo "Running comparison test, isothermal atmosphere: \n"
	@cd tests/f09comparison/                                                &&\
	../../../BART/modules/transit/transit/transit -c iso_emission.trc       &&\
	../../../BART/modules/transit/transit/transit -c iso_transmission.trc
	@echo "Isothermal comparison test complete. \n"

comparison_noinv:
	@echo "Running comparison test, noninverted atmosphere: \n"
	@cd tests/f09comparison/                                                &&\
	../../../BART/modules/transit/transit/transit -c noinv_emission.trc     &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   noinv_transmission.trc
	@echo "Noninverted comparison test complete. \n"

comparison_inv:
	@echo "Running comparison test, inverted atmosphere: \n"
	@cd tests/f09comparison/                                                &&\
	../../../BART/modules/transit/transit/transit -c inv_emission.trc       &&\
	../../../BART/modules/transit/transit/transit -c inv_transmission.trc
	@echo "Inverted comparison test complete. \n"

retrievals:
	@echo "Running retrievals: \n"

retrieval_iso:
	@echo "Running retrievals, isothermal atmosphere: \n"
	@cd tests/r01isothermal/                                                &&\
	../../../BART/BART.py -c iso_emission.brt

retrieval_noinv:
	@echo "Running retrievals, noninverted atmosphere: \n"
	@cd tests/r02noninverted/                                               &&\
	../../../BART/BART.py -c noinv_emission.brt

retrieval_inv:
	@echo "Running retrievals, inverted atmosphere: \n"
	@cd tests/r03inverted/                                                  &&\
	../../../BART/BART.py -c inv_emission.brt

retrieval_WASP12b:
	@echo "Running retrieval, WASP-12b: \n"

plots:
	@echo "Making plots..."
	@cd lib/ && makeplots.py
	@echo "Plotting complete.\n"

fin:
	@echo "Done!\n"

clean:
	@cd code-output/01BART/             &&\
	rm -f *.dat *.txt *.opt *.npz
	@cd results/01BART/                 &&\
	rm -f *.png


