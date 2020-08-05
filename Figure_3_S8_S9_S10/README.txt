This folder contains MATLAB files needed to generates Figure 3 and S8-S10. 

The files needed to run the Jupiter Notebook:  

plot_figure_3_S8_S9_S10.ipynb

are supplied here (since they take a long time to make!):

- SAMPLES_DI5_15_100000.npz
- SAMPLES_DI15_25_100000.npz
- SAMPLES_syn_DI5_15_100000.npz
- SAMPLES_syn_DI15_25_100000.npz

Running this script will produce separate files for each figure.

To regenerate each of these numpy arrays, run the following python scripts:

- pylfire_real_data.py
- pylfire_real_data2.py
- pylfire_syn_data.py
- pylfire_syn_data2.py

To run these scripts, you need to first install the software package ELFI (Engine for Likelihood-Free Inference).  Instructions can be found here:

https://elfi.readthedocs.io/en/latest/installation.html 

These python scripts call the MATLAB files:
- gas_model.m (this runs simulations in parallel)
- multistrain_model42.m (this calls the simulator once and generates summary statistics from the simulator output)
- simulator2.m (this runs the model once with an additional observation process that mimics the longitudinal data collection protocol in McDonald et al., 2008)
- parameters.m (which defines core model parameters)
- initialise_agents.m (which initialises the population in each simulation)
- pregenerate_random_numbers.m (which generates random numbers for the simulator)



