This folder contains MATLAB files needed to generate Figure 2.  

The mat files needed to run: 

generate_figure2.m 

are supplied here:

figure2_ABC.mat
figure2_DEF.mat
figure2_GHI.mat

Running this script will produce separate files for each figure panel and save them in the folder: 

figure2_figure_files

To regenerate these mat files, run: 

generate_data_for_figure2.m

This calls the following files from the Main_code_files directory:
	- parameters.m (which defines core model parameters)
	- initialise_agents.m (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator.m (which simulates transmission)
 
