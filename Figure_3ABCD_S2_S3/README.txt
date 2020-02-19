This folder contains MATLAB files needed to generate Figure 3, S2, and S3.  

The mat files needed to run: 

generate_figure3.m 

are supplied here (since they take a long time to make!):

figure3_AB.mat
figure3_CD.mat
figureS2_AB.mat
figureS2_CD.mat

Running this script will produce separate files for each figure panel and save them in the folder: 

figure3_figure_files

To regenerate these mat files, run

generate_data_for_figure3.m  

This calls the following files from the Main_code_files directory:
	- parameters.m (which defines core model parameters)
	- initialise_agents. (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator.m (which simulates transmission)

 
