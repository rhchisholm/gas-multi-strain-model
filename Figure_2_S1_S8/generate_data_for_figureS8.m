% Produces mat files that are used in generate_figure2.m to generate
% Figure S8 in Chisholm et al., Unravelling the immune response to Group A 
% Streptococcus infection from population-level observations of 
% prevalence and strain diversity

   
% Transmission scenarios: 
% Scenario 1: Figure S8 panels A-C, 
% Scenario 2: Figure S8 panels D-F,
% Scenario 3: Figure S8 panels G-I.

clear all
close all

addpath('../Main_code_files')

if ~exist('mat_files', 'dir')
    mkdir('mat_files')
end

for Scenario = 1 : 3
    
    if Scenario == 1
        sigma = 0;
        omega = 0;
    elseif Scenario == 2
        sigma = 0.5;
        omega = 0.1;
    else
        sigma = 0.5;
        omega = 0.5;
    end

    rng(2) % seed for the random numbers, needed to reproduce a simulation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Duration of simulation (years)
    DurationSimulation = 10;
    % Number of strains
    Nstrains = 42;
    % Total number of hosts
    Nagents = 2500; 
    % Age that all hosts die (years)
    AgeDeath = 71; 
    % Basic Reproduction Number 
    BasicReproductionNumber = 2.07;  
    % Duration of immunity (weeks)
    Dimmunity = 0.5 * 52.14; 
    % Resistance to co-infection
    x = 10;
    % Migration rate per week per population
    alpha = 3;  % daily per capita rate is alpha / 7 / Nagents
    % Number of contacts per week
    Cperweek = 34.53;
                   
    params = double([DurationSimulation, Nstrains, Dimmunity, ...      
                sigma, omega, x, ...              
                Cperweek, Nagents, alpha, ... 
                AgeDeath, BasicReproductionNumber]);                               
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialise agents %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [AgentCharacteristics, ImmuneStatus, ~] = initialise_agents(params);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Simulate transmission %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [SSPrev,AgentsInfectedByKStrains] = ...
                    simulator(AgentCharacteristics, ImmuneStatus, params, 0, 1);

    if Scenario == 1
        save('mat_files/figureS8_ABC.mat');
    elseif Scenario == 2
        save('mat_files/figureS8_DEF.mat');
    else
        save('mat_files/figureS8_GHI.mat');
    end
    
end



