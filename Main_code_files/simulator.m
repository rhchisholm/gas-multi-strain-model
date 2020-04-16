function [SSPrev,AgentsInfectedByKStrains] = ...
    simulator(AgentCharacteristics, ImmuneStatus, params, specifyPtransmission, cross_immunity_effect_on_coinfections)

% cross_immunity_effect_on_coinfections : 1 is on, 0 is off

[Nagents , ~, Nst, AgeDeath, ~, ~,...
    Cpertimestep, MRpertimestep, Precovery, Pimmunityloss, ...
    Ptransmission, x, StrengthImmunity, Immunity, ...
    StrengthCrossImmunity, prevalence_in_migrants, CCC,...
    ~, Ntimesteps, dt_years] = ...
        parameters(params);

pregenerate_random_numbers;

if specifyPtransmission == 1
    Ptransmission = 0.0301;
end

% Include effect that for co-infected hosts that gain immunity to a
% particular strain, the residual infections (of other strains) will clear
% faster due to cross immunity
dt= dt_years * 52.14;
Rrecovery = -log(1 - Precovery)/dt; % calculate base recovery rate from base probability of recovery
if StrengthCrossImmunity ~= 1
    Rrecovery_cici = 1 / ((1 / Rrecovery) * (1 - StrengthCrossImmunity)); % modify rate so that clearance rate speeds up due to cross immunity 
    Precovery_cici = 1 - exp(- dt * Rrecovery_cici); % calculate probability from rate
else
    Precovery_cici = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Store prevalence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, Nstrains] = size(ImmuneStatus);
% Prevalence of each strain through time
SSPrev = zeros(Nstrains, Ntimesteps);

BB = AgentCharacteristics(:,1:end-1);
SSPrev(:,1) = sum(BB,1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Store co-infection distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AgentsInfectedByKStrains = zeros(Nstrains, Ntimesteps);

if sum(sum(BB))>1
    counts = sum(BB,2);
    counts(counts==0)=[];
    [counts,K] = hist(counts,unique(counts)); % Counts number of agents (counts) with K strains
elseif sum(sum(BB))==1
    counts = 1; 
    K = 1;
else
    counts = 0; 
    K = 1;
end

AgentsInfectedByKStrains(K,1) = counts;

clear BB counts K

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MaxStrainNumber = Nst;

KeepIndexAllAgents = 1 : Nagents;

CICI = zeros(size(AgentCharacteristics(:,1:MaxStrainNumber)));

for i = 1 : Ntimesteps - 1  

    % Store current agents' characteristics here
    CurrentCharacteristics = AgentCharacteristics;
    CurrentImmune = ImmuneStatus;
    % Store current agents' characteristics without age data here
    DD = CurrentCharacteristics(:,1:MaxStrainNumber);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% RECOVERY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Find index of agents that are infected.   
    %IAgentsInf = find(DD > 0);    
    % Separate into the co-infected cross immune and otherwise:
    IAgentsInf_CICI = find(DD > 0 & CICI > 0);
    IAgentsInf_normal = find(DD > 0 & CICI == 0);

    % If n>=1 copies of an infection in host, and host clears one
    % copy, the host clears all copies.  
    % Remove index of agents from IAgentsInf that do not recover    
    %IAgentsInf(rand(size(IAgentsInf)) < (1 - Precovery)) = [];
    IAgentsInf_normal(rand(size(IAgentsInf_normal)) < (1 - Precovery)) = [];
    IAgentsInf_CICI(rand(size(IAgentsInf_CICI)) < (1 - Precovery_cici)) = [];
    
    % Update infection and immune status
    %AgentCharacteristics(IAgentsInf) = 0;
    %ImmuneStatus(IAgentsInf) = 1 * Immunity;
    AgentCharacteristics(IAgentsInf_normal) = 0;
    AgentCharacteristics(IAgentsInf_CICI) = 0;
    CICI(IAgentsInf_CICI) = 0;
    % Only normal recoveries gain strain-specific immunity
    ImmuneStatus(IAgentsInf_normal) = 1 * Immunity;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% WANING IMMUNITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Find index of agents that have immunity
    IAgentsImm = find(CurrentImmune == 1);

    % Loss of immunity occurs with probability Pimmunityloss. 
    % Remove index of agents from IAgentsImm that do not lose immunity.
    IAgentsImm(rand(size(IAgentsImm)) < (1 - Pimmunityloss)) = [];

    % Update agents that lose immunity  
    ImmuneStatus(IAgentsImm) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TRANSMISSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % Find row indices of agents that are infected. And 
    % store in RowIndexInfected.
    % Start off with all agents and remove those which are not infected
    RowIndexInfected = KeepIndexAllAgents'; % all agents
    G = sum(DD,2);
    RowIndexInfected(G==0)=[];

    if ~isempty(RowIndexInfected)  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine susceptibility of each agent to each strain:    
        % InfectionProbability stores the probability that each agent
        % will contract each strain following a contact

        % Infection status of agent:  DD
        % Number of co-infections of each agent:
        TotalInfections = sum(DD,2);

        % Immunity status of agent: CurrentImmune
        % Account for resistance to co-infection in Ptransmission1
        Ptransmission1 = Ptransmission * ...
                            (1 - TotalInfections / CCC) .^ (x);  
        Ptransmission1 = repmat(Ptransmission1,1,Nstrains);
        InfectionProbability = Ptransmission1;
        
        % Account for ss-immunity in InfectionProbability:
        if StrengthImmunity > 0
            InfectionProbability(CurrentImmune == 1) = ...
             Ptransmission1(CurrentImmune == 1) * (1 - StrengthImmunity);
        end

        % Account for cs-immunity in InfectionProbability:
        if StrengthCrossImmunity > 0
            CurrentImmunityToAnyStrain = repmat(any(CurrentImmune == 1,2),1,Nstrains);
            InfectionProbability(CurrentImmune == 0 & CurrentImmunityToAnyStrain == 1) = ...
                Ptransmission1(CurrentImmune == 0 & CurrentImmunityToAnyStrain == 1) *...
                    (1 - StrengthCrossImmunity); 
        end

         % For each agent that is infected
         for j = 1 : length(RowIndexInfected)
             
            % Here are the strains infecting the agent:
            InfectingStrainsAgent = find(DD(RowIndexInfected(j),:));

            % Determine number of contacts with other agents
            X = ContactRand(countCR,1);
            countCR = countCR + 1;
            if countCR > length(ContactRand)
                ContactRand = poissrnd(Cpertimestep,1e6,1);
                countCR = 1;
            end

            % If the agent makes contact with other agents 
            if X > 0
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Find out which agents are contacted. Here,
                % IndexOfContacts = orignal row indices of contacted agents
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Sample X agents with replacement (using 
                % pregenerated random numbers)from all agents
                IndexOfContacts = ceil(Nagents*SamplingContactsRand(countSCR:countSCR+X-1,1));
                % Update counter for pregenerated random numbers:
                countSCR = countSCR+X;
                % Regenerate random numbers if they are all used up:
                if countSCR > length(SamplingContactsRand) - Nagents
                    SamplingContactsRand = rand(1e6,1);
                    countSCR = 1;
                end
                % if this sample includes the infected agent of interest,
                % remove these and resample (the slow way):
                if ~isempty(IndexOfContacts(IndexOfContacts==RowIndexInfected(j)))
                    sampleextra = ...
                        datasample(KeepIndexAllAgents(KeepIndexAllAgents~=RowIndexInfected(j)),...
                                length(IndexOfContacts(IndexOfContacts==RowIndexInfected(j))));
                    IndexOfContacts = [IndexOfContacts(IndexOfContacts~=RowIndexInfected(j));...
                                        sampleextra'];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Determine whether transmission of any strains occurs to 
                % any of these susceptible contacts. Transmission Rule: 
                % choose strain that will transmit to each host randomly
                % then determine whether transmission occurs. More than 
                % one transmission event can occur in the time step, but 
                % susceptible hosts may only have one infection per time 
                % step 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Choose one strain to attempt transmission to each host 
                % randomly:
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % % SLOW WAY: 
                %  if length(InfectingStrainsAgent) == 1
                %      InfectingStrains = datasample(InfectingStrainsAgent,X);
                %  else
                %      InfectingStrains = datasample(InfectingStrainsAgent,X)';
                %  end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % FAST WAY:
                InfectingStrains = ceil(length(InfectingStrainsAgent)*SamplingContactsRand(countSCR:countSCR+X-1,1));
                
                if length(InfectingStrainsAgent) == 1
                    InfectingStrains = InfectingStrainsAgent(InfectingStrains);
                else
                    InfectingStrains = InfectingStrainsAgent(InfectingStrains)';
                end
                
                % Update counter for pregenerated random numbers:
                countSCR = countSCR+X;
                
                % Regenerate random numbers if they are all used up:
                if countSCR > length(SamplingContactsRand) - Nagents
                    SamplingContactsRand = rand(1e6,1);
                    countSCR = 1;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Susceptibility of contacts to infecting strains:
                %idx = sub2ind(size(InfectionProbability), IndexOfContacts,InfectingStrains); % <- SLOW
                idx = IndexOfContacts + (InfectingStrains - 1) * size(InfectionProbability,1); % <- FAST
                SusceptibilityStatusContacts =  InfectionProbability(idx); 
                
                % Determine successful transmissions:
                NewInfections = SamplingContactsRand(countSCR:countSCR+X-1,1)<SusceptibilityStatusContacts;
                
                % Update counter for pregenerated random numbers:
                countSCR = countSCR+X;
                
                % Regenerate random numbers if they are all used up:
                if countSCR > length(SamplingContactsRand) - Nagents
                    SamplingContactsRand = rand(1e6,1);
                    countSCR = 1;
                end
                
                if any(NewInfections)
                    
                    % Find new infections:
                    IndexOfContacts = IndexOfContacts(NewInfections == 1);
                    InfectingStrains = InfectingStrains(NewInfections == 1);
                    
                    % If one contact is infected more than once, remove 
                    % extra successful transmissions:
                    % % SLOW:
                    %[IndexOfContacts, temp, ~] = unique(IndexOfContacts); %IndexOfContacts now sorted, does this matter?
                    % InfectingStrains = InfectingStrains(temp);
                    % FAST:
                    [IndexOfContacts, temp] = sort(IndexOfContacts);
                    InfectingStrains = InfectingStrains(temp);
                    InfectingStrains = InfectingStrains([true;diff(IndexOfContacts(:))>0]);
                    IndexOfContacts=IndexOfContacts([true;diff(IndexOfContacts(:))>0]);
                    
                    % Update infection status:
                    %idx = sub2ind(size(AgentCharacteristics), IndexOfContacts, InfectingStrains); <- SLOW
                    idx = IndexOfContacts + (InfectingStrains - 1) * size(AgentCharacteristics,1); % <- FAST
                    AgentCharacteristics(idx) = CurrentCharacteristics(idx) + 1;
                    
                    if cross_immunity_effect_on_coinfections == 1
                        % Update matrix that tracks fast recoveries (CICI=1, fast recover, CICI=0, normal)
                        % Build tempAC2, which represents the additional infections that will 
                        % now resolve faster due to cross immunity, that need
                        % to be added to CICI.
                        tempAC1 = AgentCharacteristics(:,1:Nstrains);
                        tempAC2 = zeros(size(tempAC1));
                        % remove all infections of newly acquired strains from
                        % tempAC1
                        tempAC1(idx) = 0;
                        % For newly infected agents infected with other
                        % strains, make tempAC2=1 for those strains
                        tempAC1 = tempAC1(IndexOfContacts,:);  
                        tempAC1(tempAC1>0)=1;
                        tempAC2(IndexOfContacts,:)=tempAC1;
                        % tempAC2 now represents additional infections that will 
                        % now resolve faster, need to add this to CICI: 
                        CICI=CICI+tempAC2;
                        CICI(CICI>1)=1;
                    end
                    
                    % Make infected susceptibles not susceptible for next 
                    % infected agent of interest:
                    InfectionProbability(IndexOfContacts,:) = ...
                        zeros(length(IndexOfContacts),Nstrains);
                
                end
                
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% AGING, BIRTH & DEATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % Aging:
    AgentCharacteristics(:,Nstrains + 1) = dt_years + ...
                    CurrentCharacteristics(:,Nstrains + 1);

    % Death and Birth: remove agents that are older than 
    % AgeDeath, replace with susceptible newborn (age 0.001)
    D = find(AgentCharacteristics(:,Nstrains + 1) > AgeDeath);
    AgentCharacteristics(D,1:Nstrains) = 0;
    ImmuneStatus(D,:) = 0;
    AgentCharacteristics(D,Nstrains + 1) = 0.001;
    CICI(D,:)=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MIGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
    % Remove agents randomly, replace with agents with
    % random age, new infection state, and no immune history
    
    % Determine number of migrants from pregenerated random numbers:
    NumMig = MRRand(countMR,1);
    % Update counter for pregenerated random numbers:
    countMR = countMR + 1;
    % Regenerate random numbers if they are all used up:
    if countMR > length(MRRand)
        MRRand = poissrnd(MRpertimestep,1e6,1);
        countMR = 1;
    end
    
    % Agents leaving system:
    D = randperm(Nagents,NumMig)';
    
    % Determine characteristics of incoming migrants
    % Their infection status:
    Infected_migrants = (rand(size(D)) < prevalence_in_migrants);
    % Number of infected migrants:
    NumIMig = sum(Infected_migrants);
    
    % Determine migrant strains
    %MigrantStrains = randperm(MaxStrainNumber,NumIMig);
    %MigrantStrains = MigrantStrains(randperm(length(MigrantStrains)));
    MigrantStrains = randi(MaxStrainNumber,1,NumIMig);
        
    countM = 1; % count for immigrants into all populations
    countIM = 1; % count for infected immigrants into all populations
    % for each population:
        
    % Update infection, immune, age status of the population:
    for h = 1 : NumMig
        % new migrants have no immunity:
        ImmuneStatus(D(h),:) = zeros(1,Nstrains);
        CICI(D(h),:)=0;
        % new migrants have random age:
        AgentCharacteristics(D(h),end) = rand * AgeDeath;
        % new migrants have up to one infection:
        AgentCharacteristics(D(h),1:end-1) = zeros(1,Nstrains);  
        if Infected_migrants(countM)>0
            AgentCharacteristics(D(h),MigrantStrains(countIM)) = 1;                        
            countIM = countIM + 1;
        end

        countM = countM + 1;
            
    end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREVALENCE and COINFECTION DISTRIBUTIONS %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    BB = AgentCharacteristics(:,1:end-1);
    SSPrev(:,i + 1) = sum(BB,1)';

    if sum(sum(BB))>1
        counts = sum(BB,2);
        counts(counts==0)=[];
        [counts,K] = hist(counts,unique(counts)); % Counts number of agents (counts) with K strains
    elseif sum(sum(BB))==1
        counts = 1; 
        K = 1;
    else
        counts = 0; 
        K = 1;
    end

    AgentsInfectedByKStrains(K,i+1) = counts;
    
end

