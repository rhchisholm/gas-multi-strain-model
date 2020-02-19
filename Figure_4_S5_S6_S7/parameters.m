function [Nagents , Nstrains,  Nst, AgeDeath, NI0perstrain, NR0perstrain,...
    Cpertimestep, MRpertimestep, Precovery, Pimmunityloss, ...
    Ptransmission, x, StrengthImmunity, Immunity,...
    StrengthCrossImmunity, prevalence_in_migrants, CCC,...
    time, Ntimesteps, dt_years] = parameters(param)

multiplier = 1; %(1+1) * (param(9)>0) + 1 * (param(9)==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Duration of simulation (weeks)
Endtime = param(1) * 52.14; 
% Time step duration (weeks)
dt = 1/7; 
% Time step duration (years)
dt_years = dt / 52.14;
% Vector of time steps
time = 0:dt:Endtime;
% Total number of timesteps
Ntimesteps = length(time); 
% Initial number of infections of each strain
NI0perstrain = 10; 
% Initial number of agents immune to each strain
NR0perstrain = 10; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Epidemiological Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of initial strains  
Nstrains = param(2) * multiplier;
Nst = param(2);
% Total number of hosts
Nagents = param(8); 
% Age that all hosts die (years)
AgeDeath = param(10); 
% Duration of infection (weeks)
Dinfection = 2; 
% Duration of immunity (weeks)
Dimmunity = param(3);
% Co-infection carrying capacity
CCC =  param(2);
% Resistance to co-infection
x = param(6);
% Number of migrations per week
MR = param(9);
% Proportion of infected migrants
prevalence_in_migrants = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Strain-specific and cross-strain immunity assumptions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strength of strain-specifc immunity 
StrengthImmunity = param(4);
if StrengthImmunity>1
    StrengthImmunity = 1;
elseif StrengthImmunity<0
    StrengthImmunity = 0;
end
% Strain-specific immunity: 0 = no  s-s immunity, 1 = waning immunity 
Immunity = (StrengthImmunity > 0);
% Strength of cross immunity 
StrengthCrossImmunity = param(5);
if StrengthCrossImmunity > 1
    StrengthCrossImmunity = 1;
elseif StrengthCrossImmunity<0
    StrengthCrossImmunity = 0;
end
% Cross immunity: 0 = no cross immunity, 1 = waning immunity
%Crossimmunity = (StrengthCrossImmunity > 0); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of rates and probabilities %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rates per unit time (1/weeks)
Rrecovery = 1 / Dinfection; % recovery
Rimmunityloss = 1 / Dimmunity; % waning immunity

% Probabilities per timestep (assuming exponential distribution)
Precovery = 1 - exp(- dt * Rrecovery); % recovery
Pimmunityloss = 1 - exp(- dt * Rimmunityloss); % waning immunity

% Death / birth rate
Rdeath = 1 / AgeDeath / 52.14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Migration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mean number of hosts migrating (per timestep)
MRpertimestep = MR * dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Contacts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean number of contacts (per timestep)
Cpertimestep = param(7) * dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Reproduction Number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BasicReproductionNumber =  Ptransmission * param(7) / ...
                            %(Rdeath + Rrecovery + MR / Nagents);

% Base probability of transmission per contact
Ptransmission = (Rdeath + Rrecovery + MR / Nagents) * ...
                param(11) / param(7);

