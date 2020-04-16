This folder contains MATLAB files needed to generate Figure 2, S1 and S8.  

The mat files needed to run: 

generate_figure2.m
generate_figureS1.m
generate_figureS8.m 

are supplied here:

figure2_ABC.mat
figure2_DEF.mat
figure2_GHI.mat
figureS1_ABC.mat
figureS1_DEF.mat
figureS1_GHI.mat
figureS8_ABC.mat
figureS8_DEF.mat
figureS8_GHI.mat

Running these scripts will produce separate files for each figure panel and save them in the folder: 

figure_files

To regenerate these mat files, run: 

generate_data_for_figure2.m
generate_data_for_figureS1.m
generate_data_for_figureS8.m

These call the following files from the Main_code_files directory:
	- parameters.m (which defines core model parameters)
	- initialise_agents.m (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator.m (which simulates transmission)
 
