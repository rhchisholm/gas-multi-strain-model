This folder contains MATLAB files needed to generate Figure 2A-D, S2, S3, 
S6A-D and S7A-B.  

The mat files needed to run: 

generage_figure2ABCD_S2_S3_S6ABCD_S7AB.m 

are supplied here (since they take a long time to make!):

figure2_AB.mat
figure2_CD.mat
figureS2_AB.mat
figureS2_CD.mat
figureS6_AB.mat
figureS6_CD.mat
figureS7_AB.mat
figureS7_CD.mat

Running this script will produce separate files for each figure panel and save them in the folder: 

figure_files

To regenerate these mat files, run

generate_data_for_figure2ABCD_S2_S3_S6ABCD_S7AB.m  

This calls the following files from the Main_code_files directory:
	- parameters.m (which defines core model parameters)
	- initialise_agents. (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator.m (which simulates transmission)

 
