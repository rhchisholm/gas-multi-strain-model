% Produces mat files that are used in generate_figure3ABCD_S2_S3.m to generate
% Figure 3A-D, S2 and S3 in Chisholm et al., Unravelling the immune response to Group A 
% Streptococcus infection from population-level observations of 
% prevalence and strain diversity

addpath('../Main_code_files')

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool;
end

if ~exist('mat_files', 'dir')
    mkdir('mat_files')
end

% Scenario 1: data for Figure 3, A-B
% Scenario 2: data for Figure 3, C-D
% Scenario 3: data for Figure S2, A-B
% Scenario 4: data for Figure S2, C-D
% Scenario 5: data for Figure S9, A-B
% Scenario 6: data for Figure S9, C-D
% Scenario 7: data for Figure S10B
% Scenario 8: data for Figure S10B

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Duration of simulation (years)
DurationSimulation = 30;
% Number of time steps
Ntimesteps = length(0:1/7:DurationSimulation*52.14);
% Number of strains
Nstrains = 42;
% Total number of hosts
Nagents = 2500; 
% Age that all hosts die (years)
AgeDeath = 71; 
% Basic Reproduction Number 
BasicReproductionNumber = 2.07;  
% Migration rate per week per population
alpha = 3; % daily per capita rate is alpha / 7 / Nagents
% Number of contacts per week
Cperweek = 34.53;

Dimmunityy = [0.5 5 0.5 5 0.5 5 0.5 5];
xx = [10 10 100 100 10 10 100 100];
ciss = [0 0 0 0 1 1 1 1];

for Scenario = 1:8
    
    Dimmunity = Dimmunityy(Scenario);
    x = xx(Scenario);
    cis = ciss(Scenario);
    
    % Duration of immunity (weeks)
    DI = Dimmunity * 52.14;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Strength of strain-specific immunity
    StrainSpecific = linspace(0,1,11);
    % Strength of cross-specific immunity
    CrossStrain = linspace(0,1,11);
    % Number of simulations to average over
    NumberSimulations = 10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Store simulation data here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Final prevalence of each strain
    TimePrev = zeros(length(StrainSpecific),...
                length(CrossStrain),...
                Nstrains,...
                NumberSimulations);

    % Final number of hosts with K infections       
    TimeAgentsInfectedByKStrains = zeros(length(StrainSpecific),...
                length(CrossStrain),...
                Nstrains,... 
                NumberSimulations);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Vary parameters here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1 : length(StrainSpecific)

        SSSI = StrainSpecific(i);

        for j = 1 : i

            [Scenario i j]

            SCSI = CrossStrain(j);

            A = zeros(Nstrains,Ntimesteps,NumberSimulations);
            B = zeros(Nstrains,NumberSimulations);

            params = double([DurationSimulation, Nstrains, DI, ...      
                    SSSI, SCSI, x, ...              
                    Cperweek, Nagents, alpha, ... 
                    AgeDeath, BasicReproductionNumber]);                               

            parfor k = 1 : NumberSimulations

                [AgentCharacteristics, ImmuneStatus, ~] = ...
                        initialise_agents(params);

                [output1,output2] = ...
                        simulator(AgentCharacteristics, ImmuneStatus, params, 0, cis);

                A(:,:,k) = output1;
                B(:,k) = squeeze(output2(:,end));

            end
            
            TimePrev(i,j,:,:) = squeeze(sum(A(:,end-10*365+1:1:end,:),2));
            TimeAgentsInfectedByKStrains(i,j,:,:) = B;

        end

        if Scenario == 1
            save('mat_files/figure3_AB.mat');
        elseif Scenario == 2
            save('mat_files/figure3_CD.mat');
        elseif Scenario == 3
            save('mat_files/figureS2_AB.mat');
        elseif Scenario == 4
            save('mat_files/figureS2_CD.mat');
        elseif Scenario == 5
            save('mat_files/figureS9_AB.mat');
        elseif Scenario == 6
            save('mat_files/figureS9_CD.mat');
        elseif Scenario == 7
            save('mat_files/figureS10_AB.mat');
        elseif Scenario == 8
            save('mat_files/figureS10_CD.mat');
        end

    end

    clear output1 output2

end

