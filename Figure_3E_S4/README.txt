This folder contains MATLAB files needed to generate Figure 3E and S4.  

The mat files needed to run: 

generate_figure3E_S4A.m 

are supplied here (since they take a long time to make!):

figure3E.mat
figureS4A.mat

Running this script will produce separate files for each figure panel and save them in the folder: 

figure3E_S4_figure_files

To regenerate these mat files, run

generate_data_for_figure3E_S4A.m
	
This calls the following files from the Main_code_files directory:
	- parameters.m (which defines core model parameters)
	- initialise_agents.m (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator.m (which simulates transmission)
 
The mat file needed to run: 

generate_figureS4B.m 

is supplied here (since it takes a long time to make!):

figureS4B.mat

Running this script will produce the figure panel and save it in the folder: 

figure3E_S4_figure_files

To regenerate this mat file, run

generate_data_for_figure4SB.m
	
This calls the following files from the Main_code_files directory:
	- parameters.m (which defines core model parameters)
	- initialise_agents.m (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator_TBI.m (which simulates transmission and stores time between infection 				by the same strain)
	- run_model_once_TBI.m (which calls simulator_TBI.m for given values of duration 				of immunity and R0)
