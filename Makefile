.PHONY: all some tests_noiso test_plusiso comparison synthretrievals WASP12b bart bart_LineEtal oneline fewline multiline abundance broadening blending multicia isothermal comparison_tli comparison_iso comparison_noinv comparison_inv retrieval_iso retrieval_noinv retrieval_inv WASP12b_tli WASP12b_retrieval WASP12b_LineEtal WASP12b_StevensonEtal plots fin clean

forwardmodels: hitran_linelists oneline fewline multiline broadening abundance blending multicia isothermal comparison plots fin

tests_noiso: oneline fewline multiline broadening abundance blending multicia plots fin

tests_plusiso: oneline fewline multiline broadening abundance blending multicia isothermal plots fin

comparison: comparison_tli comparison_iso comparison_noinv comparison_inv fin

synthretrievals: comparison_tli retrieval_iso retrieval_noinv retrieval_inv fin

WASP12b: WASP12b_tli WASP12b_retrieval fin


bart:
	@echo "\nCloning BART..."
	git clone --recursive https://github.com/exosports/BART ../BART/
	@echo "Finished cloning BART to a directory parallel to BARTTest.\n"
	@echo "Compiling BART..."
	@cd ../BART/modules/transit/ && make
	@cd ../BART/modules/MCcubed/ && make
	@echo "Finished compiling BART.\n"

bart_LineEtal:
	@echo "\nCloning BART_LineEtal..."
	git clone --recursive https://github.com/exosports/BART ../BART_LineEtal/
	@echo "Finished cloning BART_Line to a directory parallel to BARTTest.\n"
	@echo "Compiling BART_LineEtal..."
	@echo "Modifying transit's makesample.c and opacity.c..."
	@cp -f tests/00inputs/transit_LineEtal/*.c                                \
	                              ../BART_LineEtal/modules/transit/transit/src/.
	@cd ../BART_LineEtal/modules/transit/ && make
	@cd ../BART_LineEtal/modules/MCcubed/ && make
	@echo "Finished compiling BART_LineEtal.\n"

hitran_linelists:
	@echo "Downloading HITRAN line lists...\n"
	@cd tests/00inputs/par/                                                 &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITEMP-CO2.txt    &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITEMP-CO.txt     &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITEMP-H2O.txt    &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITRAN-CH4.txt    &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITRAN-HCN.txt    &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITRAN-C2H2.txt   &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITRAN-C2H6.txt   &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITRAN-H2.txt     &&\
	wget --user=HITRAN --password=getdata -N -i wget-list_HITRAN-NH3.txt
	@echo "Extracting archives...\n"
	@cd tests/00inputs/par/                                                 &&\
	unzip '01_*HITEMP2010.zip'                                              &&\
	unzip '02_*HITEMP2010.zip'                                              &&\
	unzip '05_*HITEMP2010.zip'                                              &&\
	unzip '06_hit12.zip'                                                    &&\
	unzip '11_hit12.zip'                                                    &&\
	unzip '23_hit08.zip'                                                    &&\
	unzip '26_hit12.zip'                                                    &&\
	unzip '27_hit12.zip'                                                    &&\
	unzip '45_hit12.zip'                                                    &&\
	rm -f *.zip
	@echo "Finished retrieving HITRAN line lists.\n"

oneline:
	@mkdir -p "code-output/01BART/f01oneline"
	@echo "Running oneline test...\n"
	@cd tests/f01oneline/                                                   &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   oneline.plc                                                          &&\
	../../../BART/modules/transit/transit/transit -c oneline_emission.trc
	@echo "oneline test complete.\n"

fewline:
	@mkdir -p "code-output/01BART/f02fewline"
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
	@mkdir -p "code-output/01BART/f03multiline"
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
	@mkdir -p "code-output/01BART/f04broadening"
	@echo "Running broadening test...\n"
	@cd tests/f04broadening/                                                &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   broadening.plc                                                       &&\
	../../../BART/modules/transit/transit/transit -c broadening_emission.trc
	@cd lib/ && python -c 'import voigtcomp;voigtcomp.comp()'
	@echo "broadening test complete.\n"

abundance:
	@mkdir -p "code-output/01BART/f05abundance"
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
	@echo "abundance test complete.\n"

blending:
	@mkdir -p "code-output/01BART/f06blending"
	@echo "Running blending test...\n"
	@cd tests/f06blending/                                                  &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   blending.plc                                                         &&\
	../../../BART/modules/transit/transit/transit -c blending_emission.trc
	@echo "blending test complete.\n"

multicia:
	@mkdir -p "code-output/01BART/f07multicia"
	@echo "Running multicia test...\n"
	@cd tests/f07multicia/                                                  &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   multicia.plc                                                         &&\
	../../../BART/modules/transit/transit/transit -c noCIA_emission.trc     &&\
	../../../BART/modules/transit/transit/transit -c oneCIA_emission.trc    &&\
	../../../BART/modules/transit/transit/transit -c twoCIA_emission.trc
	@echo "multicia test complete.\n"

generate_tli:
	@echo "Generating TLI file from HITRAN linelists...\n"
	@cd tests/00inputs/TLI/                                                 &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c pyline.plc
	@echo "TLI generation complete. \n"

isothermal:
	@mkdir -p "code-output/01BART/f08isothermal"
	@echo "Running isothermal test...\n"
	@cd tests/f08isothermal/                                                &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   isothermal.plc                                                       &&\
	../../../BART/modules/transit/transit/transit -c isothermal_emission.trc
	@echo "isothermal test complete. \n"

comparison_tli:
	@echo "Generating TLI file for comparison/synthetic retrieval tests...\n"
	@cd tests/f09comparison                                                 &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   comparison.plc
	@echo "TLI file generated for comparison tests.\n"

comparison_iso:
	@mkdir -p "code-output/01BART/f09comparison"
	@echo "Running comparison test, isothermal atmosphere: \n"
	@cd tests/f09comparison/                                                &&\
	../../../BART/modules/transit/transit/transit -c iso_emission.trc       &&\
	../../../BART/modules/transit/transit/transit -c iso_transmission.trc
	@echo "Isothermal comparison test complete. \n"

comparison_noinv:
	@mkdir -p "code-output/01BART/f09comparison"
	@echo "Running comparison test, noninverted atmosphere: \n"
	@cd tests/f09comparison/                                                &&\
	../../../BART/modules/transit/transit/transit -c noinv_emission.trc     &&\
	../../../BART/modules/transit/transit/transit -c                          \
	   noinv_transmission.trc
	@echo "Noninverted comparison test complete. \n"

comparison_inv:
	@mkdir -p "code-output/01BART/f09comparison"
	@echo "Running comparison test, inverted atmosphere: \n"
	@cd tests/f09comparison/                                                &&\
	../../../BART/modules/transit/transit/transit -c inv_emission.trc       &&\
	../../../BART/modules/transit/transit/transit -c inv_transmission.trc
	@echo "Inverted comparison test complete. \n"

retrieval_iso:
	@echo "Running retrievals, isothermal atmosphere: \n"
	@cd tests/r01isothermal/                                                &&\
	../../../BART/BART.py -c iso_emission.brt                               &&\
	../../../BART/BART.py -c iso_transmission.brt
	@echo "Isothermal retrieval test complete. \n"

retrieval_noinv:
	@echo "Running retrievals, noninverted atmosphere: \n"
	@cd tests/r02noninverted/                                               &&\
	../../../BART/BART.py -c noinv_emission.brt                             &&\
	../../../BART/BART.py -c noinv_transmission.brt
	@echo "Noninverted retrieval test complete. \n"

retrieval_inv:
	@echo "Running retrievals, inverted atmosphere: \n"
	@cd tests/r03inverted/                                                  &&\
	../../../BART/BART.py -c inv_emission.brt                               &&\
	../../../BART/BART.py -c inv_transmission.brt
	@echo "Inverted retrieval test complete. \n"

WASP12b_tli:
	@echo "Generating TLI file for WASP-12b retrieval...\n"
	@cd tests/r04WASP12b/                                                   &&\
	../../../BART/modules/transit/pylineread/src/pylineread.py -c             \
	   WASP12b.plc
	@echo "TLI file generated for WASP-12b retrieval.\n"

WASP12b_retrieval:
	@echo "Running retrieval, WASP-12b: \n"
	@cd tests/r04WASP12b/                                                   &&\
	../../../BART/BART.py -c WASP12b.brt
	@echo "WASP-12b retrieval complete.\n"

WASP12b_LineEtal:
	@echo "Running retrieval, WASP-12b, Line et al. 2014 case: \n"
	@cd tests/r04WASP12b/                                                   &&\
	../../../BART_LineEtal/BART.py -c WASP12b_LineEtal2014.brt
	@echo "WASP-12b Line et al. 2014 retrieval complete.\n"

WASP12b_StevensonEtal:
	@echo "Running retrieval, WASP-12b, Stevenson et al. 2014 case: \n"
	@cd tests/r04WASP12b/                                                   &&\
	../../../BART/BART.py -c WASP12b_StevensonEtal2014.brt
	@echo "WASP-12b Stevenson et al. 2014 retrieval complete.\n"

plots:
	@echo "Making plots..."
	@cd lib/ && ./makeplots.py
	@echo "Plotting complete.\n"

fin:
	@echo "Done!\n"

clean:
	@cd code-output/01BART/             &&\
	rm -f *
	@cd results/01BART/                 &&\
	rm -f *


