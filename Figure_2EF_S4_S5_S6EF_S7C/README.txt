This folder contains MATLAB files needed to generate Figure 2E,F, S4A,B,C, S5, S6E,F, S7C. 

The mat files needed to run: 

generate_figure2EF_S4AB_S6EF.m
generate_figureS5.m
generate_figureS4C_S7C.m

are supplied in the folder mat_files (since they take a long time to make!)

Running these scripts will produce separate files for each figure panel and 
save them in the folder: 

figure_files

To regenerate these mat files, run

generate_data_for_figure2E_S4A_S6E.m and
generate_data_for_figure2F_S4B_S6F.m and
generate_data_for_figureS5AC.m
generate_data_for_figureS5BD.m
generate_data_for_figureS4C_S7C.m
	
These call the following files from the Main_code_files directory related to our model:
	- parameters.m (which defines core model parameters)
	- initialise_agents.m (which initialises the population)
	- pregenerate_random_numbers.m (which generates random numbers for the simulator)
	- simulator.m (which simulates transmission)
    	- simulator_TBI.m (which simulates transmission and stores time between infection) 				by the same strain)
	- run_model_once_TBI.m (which calls simulator_TBI.m for given values of duration)				of immunity and R0)

and the equivalent files for the SIRIR model:
	- parameters_base_SIRIR.m (which defines core model parameters)
	- initialise_agents_SIRIR.m (which initialises the population)
	- simulator_SIRIR.m (which simulates transmission)


 
