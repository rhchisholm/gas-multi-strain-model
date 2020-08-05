function parameters = parameters_base_SIRIR(parameters,endtime)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Base model parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model setup parameters 
parameters.T = endtime * 52.14; % Duration of simulation (weeks)
parameters.dt = 1/7; % Time step duration (weeks)
parameters.dt_years = parameters.dt / 52.14; % Time step duration (years)
parameters.time = 0 : parameters.dt : parameters.T; % vector of time steps
parameters.Ntimesteps = length(parameters.time); % Number of time steps
parameters.NI0perstrain = 10; % Initial number of infections of each strain in each population
parameters.NR0perstrain = 10;

% Population parameters
parameters.Nstrains = 42;
parameters.MaxStrains = 42; % Maximum number of strains in system
%parameters.NcontactsPerWeek = 32.6 * 7;
%parameters.Ncontacts = parameters.NcontactsPerWeek * parameters.dt; % per time step
%parameters.PopSize = 2500;  %_SIRI* ones(size(parameters.Nstrains));
%parameters.MaxPopSize = 2500; %* ones(size(parameters.Nstrains)); % Population capacity of each population
parameters.NPopulations = 1; % Number of populations
parameters.AgeDeath = 71; % years
parameters.AgeCategories = [0 5 11 18 45 65 75];

% Migration
parameters.ImmigrationRatePerCapita = 1.71 * 10^(-4); % per time step, into and out of whole system
parameters.prevalence_in_migrants = 0.1; % prevalence of infection in immigrants
%parameters.ImmigrationRatePerCapitaPerWeek = parameters.ImmigrationRatePerCapita / parameters.dt; % per week, needed for calculation of R0

% Intervention parameters
parameters.TimeIntervention = 100; % years
% vschoice = 1: vaccine strains are most prevalent
% vschoise = 0: vaccine strains are random (strains 1:numvs)
parameters.vschoice = 1; 
% catchup = 1: catch-up campaign at start of intervention
% catchup = 0: no catch-up campaign 
parameters.catchup = 0;
parameters.CampaignVaccinationCoverage = 0.9;
parameters.RoutineVaccinationCoverage = 0.9;
parameters.NumberOfVaccineStrains = 0;
parameters.VaccineStrains = 1:parameters.NumberOfVaccineStrains;%randperm(parameters.MaxStrains,parameters.NumberOfVaccineStrains);
x = 1:parameters.MaxStrains;
parameters.NVaccineStrains = x(~ismember(x,parameters.VaccineStrains));
parameters.AgeOfVacccination = 0.5; %(years)
parameters.iIntervention = floor(parameters.TimeIntervention * 52.14 * 7);

% Within-host parameters
parameters.DurationInfection = 2; % Duration of infection (weeks)
parameters.SigmaShort = 0.9; % Strength of strain-specific immunity, when immunity is temporary
parameters.OmegaShort = 0.1; % Strength of cross-strain immunity, when immunity is temporary
parameters.SigmaLong = 1; % Strength of strain-specific immunity, when immunity is life-long
parameters.OmegaLong = 0.1; % Strength of cross-strain immunity, when immunity is life-long
parameters.NumberRepeatInfections = 2; % Number of timely infections required to obtain life-long immunity

% Parameters governing rates and probabilities of events
parameters.Rrecovery = 1 / parameters.DurationInfection; % base recovery rate per week
parameters.Precovery = 1 - exp(- parameters.dt * parameters.Rrecovery); % base recovery probability per time step

% Recovery rate per week when there is some immunity
parameters.RrecoveryImmuneSSS = 1 / (parameters.DurationInfection * (1 - parameters.SigmaShort)); % recovery per week if agent has short-term strain specific immunity to that strain
parameters.RrecoveryImmuneSSL = 1 / (parameters.DurationInfection * (1 - parameters.SigmaLong)); % recovery per week if agent has life-long strain specific immunity to that strain
parameters.RrecoveryImmuneCSS = 1 / (parameters.DurationInfection * (1 - parameters.OmegaShort)); % recovery per week if agent has short-term cross-strain immunity to that strain
parameters.RrecoveryImmuneCSL = 1 / (parameters.DurationInfection * (1 - parameters.OmegaLong)); % recovery per week if agent has life-long cross-strain immunity to that strain

% Recovery probability per time step when there is some immunity
parameters.PrecoveryImmuneSSS = 1 - exp(- parameters.dt * parameters.RrecoveryImmuneSSS); % probabiltiy per time step if agent has short-term strain specific immunity to that strain
parameters.PrecoveryImmuneSSL = 1 - exp(- parameters.dt * parameters.RrecoveryImmuneSSL); % probabiltiy per time step if agent has life-long strain specific immunity to that strain
parameters.PrecoveryImmuneCSS = 1 - exp(- parameters.dt * parameters.RrecoveryImmuneCSS); % probabiltiy per time step if agent has short-term cross-strain immunity to that strain
parameters.PrecoveryImmuneCSL = 1 - exp(- parameters.dt * parameters.RrecoveryImmuneCSL); % probabiltiy per time step if agent has life-long cross-strain immunity to that strain

parameters.Rdeath = 1 / parameters.AgeDeath / 52.14; % death rate per week

parameters.CCC = 42;
 