This folder contains MATLAB files needed to generate Figure 3A-D, S2, S3, 
S9A-D and S10A-B.  

The mat files needed to run: 

generage_figure3ABCD_S2_S3_S9ABCD_S10AB.m 

are supplied here (since they take a long time to make!):

figure3_AB.mat
figure3_CD.mat
figureS2_AB.mat
figureS2_CD.mat
figureS9_AB.mat
figureS9_CD.mat
figureS10_AB.mat
figureS10_CD.mat

Running this script will produce separate files for each figure panel and save them in the folder: 

figure_files

To regenerate these mat files, run

generate_data_for_figure3ABCD_S2_S3_S9ABCD_S10AB.m  

This calls the following files from the Main_code_files directory:
	- parameters.m (which defines core model parameters)
	- initialise_agents. (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator.m (which simulates transmission)

 
