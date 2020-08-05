This folder contains MATLAB files needed to generate Figure S1.  

The mat files needed to run: 

generate_figureS1.m

are supplied here:

figureS1_ABCD.mat
figureS1_EFGH.mat
figureS1_IJKL.mat
figureS1_MNOP.mat
figureS1_QRST.mat
figureS1_UVWX.mat

Running these scripts will produce separate files for each figure panel and save them in the folder: 

figure_files

To regenerate these mat files, run: 

generate_data_for_figureS1A_L.m
generate_data_for_figureS1M_X.m

These call the following files from the Main_code_files directory:
	- parameters.m (which defines core model parameters)
	- initialise_agents.m (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator.m (which simulates transmission)
 
