% Produces mat files that are used in generate_figure3EF_S4AB_S9EF.m to 
% generate Figures 3F,S4B and S9F in Chisholm et al., Unravelling the 
% immune response to Group A Streptococcus infection from population-level 
% observations of prevalence and strain diversity

clear all

rng(1)

addpath('../Main_code_files')

if ~exist('mat_files', 'dir')
    mkdir('mat_files')
end

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
N = [200 1000 2000 3500 5000 10000];
% Age that all hosts die (years)
AgeDeath = 71; 
% Basic Reproduction Number 
BasicReproductionNumber = 2.07;  
% Migration rate per capita per week per population
delta = 3/2500;
% Number of contacts (per week)
Contacts = 34.53;
% Number of simulations
NumberSimulations = 10;

xall = [10 100 10];
cross_immunity_effect_on_coinfections = [0 0 1];
sigma = 1;
omega = 0.1;

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool(10);
end

% rtc 1: data for Figure 3F
% rtc 2: data for Figure S4B
% rtc 3: data for Figure S9F

for rtc = 1:3
    
    x = xall(rtc);
    cis = cross_immunity_effect_on_coinfections(rtc);
    
    % Scenario 1:  no immunity
    % Scenario 2:  6 months immunity
    % Scenario 3:  5 years immunity
    % Scenario 4:  lifelong immunity

    for Scenario = 1:4

        if Scenario == 1
            DI = 0 * 52.14;
        elseif Scenario == 2
            DI = 0.5 * 52.14;
        elseif Scenario ==3
            DI = 5 * 52.14;
        else
            DI = 71 * 52.14;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Store simulation data here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TimePrev = zeros(length(delta),...
                    length(N),...
                    NumberInitialStrains,...
                    NumberSimulations);

        TimeAgentsInfectedByKStrains = zeros(length(delta),...
                    length(N),...
                    NumberInitialStrains,...
                    NumberSimulations);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Vary parameters here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for i = 1 : length(delta)

            for j = 1 : length(N)

                A = zeros(NumberInitialStrains,Ntimesteps,NumberSimulations);
                B = zeros(NumberInitialStrains,NumberSimulations);

                PopSize = N(j);
                mr = delta(i) * N(j);

                params = double([DurationSimulation, NumberInitialStrains, DI, ...      
                        sigma, omega, x, ...              
                        Contacts, PopSize, mr, ... 
                        AgeDeath, BasicReproductionNumber]);

                [Scenario i j]

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
                
                if rtc == 1
                
                    if Scenario == 1
                        save('mat_files/figure3F_1.mat');
                    elseif Scenario == 2
                        save('mat_files/figure3F_2.mat');
                    elseif Scenario ==3
                        save('mat_files/figure3F_3.mat');
                    else
                        save('mat_files/figure3F_4.mat');
                    end
                    
                elseif rtc == 2
                    
                    if Scenario == 1
                        save('mat_files/figureS4B_1.mat');
                    elseif Scenario == 2
                        save('mat_files/figureS4B_2.mat');
                    elseif Scenario ==3
                        save('mat_files/figureS4B_3.mat');
                    else
                        save('mat_files/figureS4B_4.mat');
                    end
                    
                else
                    
                    if Scenario == 1
                        save('mat_files/figureS9F_1.mat');
                    elseif Scenario == 2
                        save('mat_files/figureS9F_2.mat');
                    elseif Scenario ==3
                        save('mat_files/figureS9F_3.mat');
                    else
                        save('mat_files/figureS9F_4.mat');
                    end
                                    
                end

            end

        end

        clear output1 output2 A B

        if rtc == 1
                
            if Scenario == 1
                save('mat_files/figure3F_1.mat');
            elseif Scenario == 2
                save('mat_files/figure3F_2.mat');
            elseif Scenario ==3
                save('mat_files/figure3F_3.mat');
            else
                save('mat_files/figure3F_4.mat');
            end

        elseif rtc == 2

            if Scenario == 1
                save('mat_files/figureS4B_1.mat');
            elseif Scenario == 2
                save('mat_files/figureS4B_2.mat');
            elseif Scenario ==3
                save('mat_files/figureS4B_3.mat');
            else
                save('mat_files/figureS4B_4.mat');
            end

        else

            if Scenario == 1
                save('mat_files/figureS9F_1.mat');
            elseif Scenario == 2
                save('mat_files/figureS9F_2.mat');
            elseif Scenario ==3
                save('mat_files/figureS9F_3.mat');
            else
                save('mat_files/figureS9F_4.mat');
            end

        end

    end
    
end



