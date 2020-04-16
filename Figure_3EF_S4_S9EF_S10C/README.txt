This folder contains MATLAB files needed to generate Figure 3E,F, S4A,B,C,  
S9E,F and S10C. 

The mat files needed to run: 

generate_figure3EF_S4AB_S9EF.m and
generate_figureS4C_S10C.m

are supplied in the folder mat_files (since they take a long time to make!)

Running these scripts will produce separate files for each figure panel and 
save them in the folder: 

figure_files

To regenerate these mat files, run

generate_data_for_figure3E_S4A_S9E.m and
generate_data_for_figureS4C_S10C.m
	
These call the following files from the Main_code_files directory:
	- parameters.m (which defines core model parameters)
	- initialise_agents.m (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator.m (which simulates transmission)
    - simulator_TBI.m (which simulates transmission and stores time between infection) 				by the same strain)
	- run_model_once_TBI.m (which calls simulator_TBI.m for given values of duration)				of immunity and R0)

 
