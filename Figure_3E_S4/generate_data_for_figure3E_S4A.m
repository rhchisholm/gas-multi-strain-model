% Produces mat files that are used in generate_figure3E_S4A.m to generate
% Figure 3E and S4A in Chisholm et al., Unravelling the within-host 
% dynamics of Group A Streptococcus from population-level observations of 
% prevalence and strain diversity

clear all

addpath('../Main_code_files')

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool;
end

if ~exist('figure3E_S4_mat_files', 'dir')
    mkdir('figure3E_S4_mat_files')
end

% Scenario 1: data for Figure 4
% Scenario 2: data for Figure 4SA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Duration of simulation (years)
DurationSimulation = 30;
% Number of time steps
Ntimesteps = length(0:1/7:DurationSimulation*52.14);
% Number of strains
NumberInitialStrains = 42;
% Total number of hosts
Nagents = 2500; 
% Age that all hosts die (years)
AgeDeath = 71;  
% Migration rate per week per population
alpha = 3; % daily per capita rate is alpha / 7 / Nagents
% Number of contacts (per week). 
Contacts = 16.70 * (1:0.5:10.5); %This equates to R0 bw 1:0.5:10.5
% Number of simulations per parameter combination
NumberSimulations = 10;
% Duration of immunity (weeks)
Dimmunity = [0.01 0.5*52.14 5*52.14 71*52.14];
% Strengths of immunity
sigma = 1;
omega = 0.1;

for Scenario = 1:2
    
    if Scenario == 1       
        x = 10;
    else
        x = 100;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Store simulation data here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    TimePrev = zeros(length(Dimmunity),...
                length(Contacts),...
                NumberInitialStrains,...
                NumberSimulations);

    TimeAgentsInfectedByKStrains = zeros(length(Dimmunity),...
                length(Contacts),...
                NumberInitialStrains,...
                NumberSimulations);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Vary parameters here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1 : length(Dimmunity)

        DI = Dimmunity(i);

        for j = 1 : length(Contacts)
            
            [Scenario,i,j]
            
            A = zeros(NumberInitialStrains,Ntimesteps,NumberSimulations);
            B = zeros(NumberInitialStrains,NumberSimulations);

            Cperweek = Contacts(j);
            
            params = double([DurationSimulation, NumberInitialStrains, DI, ...      
                    sigma, omega, x, ...              
                    Cperweek, Nagents, alpha, ... 
                    AgeDeath, BasicReproductionNumber]);

            parfor k = 1 : NumberSimulations

                [AgentCharacteristics, ImmuneStatus, ~] = ...
                        initialise_agents(params);

                [output1,output2] = ...
                        simulator(AgentCharacteristics, ImmuneStatus, params, 1);
                
                A(:,:,k) = output1;
                B(:,k) = squeeze(output2(:,end));

            end

            TimePrev(i,j,:,:) = squeeze(sum(A(:,end-10*365+1:1:end,:),2));
            TimeAgentsInfectedByKStrains(i,j,:,:) = B;
            
            if Scenario == 1
                save('figure3E_S4_mat_files/figure3E.mat');
            else
                save('figure3E_S4_mat_files/figureS4A.mat');
            end

        end

    end

    clear output1 output2 A B

    if Scenario == 1
        save('figure3E_S4_mat_files/figure3E.mat');
    else
        save('figure3E_S4_mat_files/figureS4A.mat');
    end

end


