% Produces mat files that are used in generate_figureS5.m to 
% generate SIRIS model data for Figure S5A,C 

clear all

addpath('Main_code_files')

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool(10);
end

if ~exist('mat_files', 'dir')
    mkdir('mat_files')
end

% Scenario 1: data for Figure S5A
% Scenario 2: data for Figure S5C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Duration of simulation (years)
DurationSimulation = 30;
% Number of time steps
Ntimesteps = length(0:1/7:DurationSimulation*52.14);
% Number of contacts (per week). 
Contacts = 16.70 * (1:0.5:10.5); %This equates to R0 bw 1:0.5:10.5
ContactsPerDay = Contacts / 7;
% Number of simulations per parameter combination
NumberSimulations = 10;
% Duration of immunity window (weeks)
Dimmunity = [3 19 104];

parameters.PopSize = 2500;  
parameters.MaxPopSize = 2500;
parameters.Ptransmission = 0.0301;

parameters = parameters_base_SIRIR(parameters,DurationSimulation);

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
                parameters.MaxStrains,...
                NumberSimulations);

    TimeAgentsInfectedByKStrains = zeros(length(Dimmunity),...
                length(Contacts),...
                parameters.MaxStrains,...
                NumberSimulations);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Vary parameters here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1 : length(Dimmunity)

        for j = 1 : length(Contacts)
            
            [Scenario,i,j]
            
            A = zeros(parameters.MaxStrains,Ntimesteps,NumberSimulations);
            B = zeros(parameters.MaxStrains,NumberSimulations);

            parameters.Ncontacts = Contacts(j);
            parameters.TimeWindow = Dimmunity(i);
            parameters.coinfection = x;

            parfor k = 1 : NumberSimulations

                AgentCharacteristics =  initialise_agents_SIRIR(parameters);

                [output1,output2] = ...
                        simulator_SIRIR(AgentCharacteristics, parameters, 0);
                
                A(:,:,k) = output1;
                B(:,k) = squeeze(output2(:,end));

            end

            TimePrev(i,j,:,:) = squeeze(sum(A(:,end-10*365+1:1:end,:),2));
            TimeAgentsInfectedByKStrains(i,j,:,:) = B;
            
            if Scenario == 1
                save('mat_files/figureS5A.mat');
            else
                save('mat_files/figureS5C.mat');
            end

        end

    end

    clear output1 output2 A B

    if Scenario == 1
        save('mat_files/figureS5A.mat');
    else
        save('mat_files/figureS5C.mat');
    end

end


