This folder contains MATLAB files needed to generate Figure S1.  

The mat files needed to run: 

generate_figureS1.m 

are supplied here:

figureS1_ABC.mat
figureS1_DEF.mat
figureS1_GHI.mat

Running this script will produce separate files for each figure panel and save them in the folder: 

figureS1_figure_files

To regenerate these mat files, run: 

generate_data_for_figureS1.m

This calls the following files from the Main_code_files directory:
	- parameters.m (which defines core model parameters)
	- initialise_agents.m (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator.m (which simulates transmission)



 
