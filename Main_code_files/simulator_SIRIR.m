function [SSPrev,AgentsInfectedByKStrains] = ...
    simulator_SIRIR(AgentCharacteristics, parameters, intervention)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Pregenerate random numbers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Generate random numbers for contact events
    ContactRand = poissrnd(parameters.Ncontacts,[1e6,1]);
    countCR = 1;

    % Generate random numbers for immigration events
    MRRand = poissrnd(parameters.ImmigrationRatePerCapita * ...
                                parameters.PopSize,1e6,1);        
    countMR = 1;

    % Generate random numbers for sampling contacts
    SamplingContactsRand = rand(1e6,1);
    countSCR = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Store Summary Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SummaryStatistics.SingleStrainPrevTime = zeros(parameters.MaxStrains,parameters.Ntimesteps);    
    SummaryStatistics.SingleStrainLongImmunityTime = zeros(parameters.MaxStrains,parameters.Ntimesteps); 
    SummaryStatistics.SingleStrainPrevTime(:,1) = sum(AgentCharacteristics.Infections,1)';
    SummaryStatistics.SingleStrainLongImmunityTime(:,1) = sum(AgentCharacteristics.ImmunityLong,1)';
    SummaryStatistics.CoinfectionDistributionTime = zeros(parameters.MaxStrains,parameters.Ntimesteps);   
    SummaryStatistics.LongImmunityDistributionTime = zeros(parameters.MaxStrains,parameters.Ntimesteps);
    SummaryStatistics.LongImmunityDistributionTimeVS = zeros(parameters.NumberOfVaccineStrains,parameters.Ntimesteps);
    SummaryStatistics.LongImmunityDistributionTimeNVS = zeros(parameters.MaxStrains - parameters.NumberOfVaccineStrains,parameters.Ntimesteps);
    SummaryStatistics.AvgVaccinationCoverage = zeros(1,parameters.Ntimesteps);

    % Calculates the number of agents (counts) with K infections
    [counts,K] = calculate_coinfection_distribution(AgentCharacteristics.Infections);  
    SummaryStatistics.CoinfectionDistributionTime(K,1) = counts;
    % Calculates the number of agents (counts) with long immunity to K
    % strains

    [counts,K] = calculate_coinfection_distribution(AgentCharacteristics.ImmunityLong);  
    SummaryStatistics.LongImmunityDistributionTime(K,1) = counts;
    [counts,K] = calculate_coinfection_distribution(AgentCharacteristics.ImmunityLong(:,parameters.VaccineStrains));  
    SummaryStatistics.LongImmunityDistributionTimeVS(K,1) = counts;
    [counts,K] = calculate_coinfection_distribution(AgentCharacteristics.ImmunityLong(:,parameters.NVaccineStrains));  
    SummaryStatistics.LongImmunityDistributionTimeNVS(K,1) = counts;

    SummaryStatistics.PopulationSizeTime = zeros(1,parameters.Ntimesteps);
    SummaryStatistics.PopulationSizeTime(1,1) = parameters.PopSize;

    SummaryStatistics.AvgLongImmunity = zeros(1,parameters.Ntimesteps);
    SummaryStatistics.AvgLongImmunity(1,1) =  calculate_avg_immunity(SummaryStatistics.LongImmunityDistributionTime(:,1)) ./ parameters.PopSize;

    SummaryStatistics.AvgVaccinationCoverage(1,1) =  mean(AgentCharacteristics.VaccinationStatus);

    SummaryStatistics.AvgLongImmunityVS = zeros(1,parameters.Ntimesteps);
    SummaryStatistics.AvgLongImmunityNV = zeros(1,parameters.Ntimesteps);
    SummaryStatistics.AvgLongImmunityVS(1,1) =  calculate_avg_immunity(SummaryStatistics.LongImmunityDistributionTimeVS(:,1)) ./ parameters.PopSize;
    SummaryStatistics.AvgLongImmunityNVS(1,1) =  calculate_avg_immunity(SummaryStatistics.LongImmunityDistributionTimeNVS(:,1)) ./ parameters.PopSize;

    SummaryStatistics.PrevalenceVaccinated = zeros(1,parameters.Ntimesteps);
    SummaryStatistics.PrevalenceUnvaccinated = zeros(1,parameters.Ntimesteps);
    SummaryStatistics.NumberVaccinated = zeros(1,parameters.Ntimesteps);
    SummaryStatistics.NumberUnvaccinated = zeros(1,parameters.Ntimesteps);

    [SummaryStatistics.PrevalenceVaccinated(1,1),...
    SummaryStatistics.PrevalenceUnvaccinated(1,1),...
    SummaryStatistics.NumberVaccinated(1,1),...
    SummaryStatistics.NumberUnvaccinated(1,1)] = ...
            calculate_prev_vacc(AgentCharacteristics.Infections,...
                                AgentCharacteristics.VaccinationStatus);

    SummaryStatistics.NumberVaccineStrainsPresent = zeros(1,parameters.Ntimesteps);
    SummaryStatistics.NumberNonVaccineStrainsPresent = zeros(1,parameters.Ntimesteps);

    [SummaryStatistics.NumberVaccineStrainsPresent(1,1),...
    SummaryStatistics.NumberNonVaccineStrainsPresent(1,1)] = ...
            calculate_number_strain_present(AgentCharacteristics.Infections,...
            parameters.VaccineStrains, parameters.NVaccineStrains);

    SummaryStatistics.AgeDistributionPrevalence = zeros(length(parameters.AgeCategories)-1,parameters.Ntimesteps);
    SummaryStatistics.AgeDistributionPrevalence(:,1) = ...
            calculate_age_dist_prev(AgentCharacteristics.Infections,AgentCharacteristics.Age,parameters.AgeCategories);   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Preallocate memory for cell arrays %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    KeepIndexAllAgents = 1 : parameters.PopSize;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1 : parameters.Ntimesteps - 1 % for each time step

        CurrentTime = i * parameters.dt; % this is the current time in weeks
        
        if i == parameters.iIntervention

            if parameters.vschoice == 1
                [~,SortedStrains] = sort(SummaryStatistics.SingleStrainPrevTime(:,i),'descend');
                parameters.VaccineStrains = SortedStrains(1:parameters.NumberOfVaccineStrains);
                x = 1:parameters.MaxStrains;
                parameters.NVaccineStrains = x(~ismember(x,parameters.VaccineStrains));
            end

            if parameters.catchup == 0 % No catch-up campaign:
                TreatedAgents = [];
            else % With catch-up campaign:
                TreatedAgents = find(AgentCharacteristics.Age > 5 & ...
                   AgentCharacteristics.Age <= 11);
                TreatedAgents(rand(size(TreatedAgents))>parameters.CampaignVaccinationCoverage)=[];
                %TreatedAgents = randperm(parameters.PopSize,...
                               % round(parameters.PopSize * parameters.CampaignVaccinationCoverage));
            end

            AgentCharacteristics.VaccinationStatus(TreatedAgents,1) = ones(length(TreatedAgents),1);
            AgentCharacteristics.ImmunityLong(TreatedAgents,parameters.VaccineStrains) = ...
                ones(length(TreatedAgents),parameters.NumberOfVaccineStrains);
            AgentCharacteristics.ImmunityShort(TreatedAgents,parameters.VaccineStrains) = ...
                zeros(length(TreatedAgents),parameters.NumberOfVaccineStrains);
            AgentCharacteristics.Infections(TreatedAgents,parameters.VaccineStrains) = ...
                zeros(length(TreatedAgents),parameters.NumberOfVaccineStrains);

        end

        % Store current agents' characteristics here
        CurrentInfections = AgentCharacteristics.Infections;
        CurrentImmuneLong = AgentCharacteristics.ImmunityLong;
        CurrentImmuneShort = AgentCharacteristics.ImmunityShort;
        CurrentRecoveryTime = AgentCharacteristics.RecoveryTime;
        CurrentVaccineCoverage = AgentCharacteristics.VaccinationStatus;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% WANING IMMUNITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Agents lose short immunity on average after parameters.TimeWindow
        % time units.

        % Find index of agents that have short term immunity and time 
        % since last recovery exceeds parameters.TimeWindow (so that 
        % they now lose short-term immunity)
        IAgentsImm = find(CurrentImmuneShort > 0 & ...
            (CurrentTime - CurrentRecoveryTime > parameters.TimeWindow));

        % Update agents that lose short-term immunity  
        AgentCharacteristics.ImmunityShort(IAgentsImm) = 0;
        AgentCharacteristics.RecoveryTime(IAgentsImm) = 0;

        % During a time-step, we assume that if an agent loses 
        % short term immunity to a particular strain, then they  
        % cannot then obtain long-term immunity to that strain 
        % during the same time  step.  Therefore, need to up update
        % CurrentImmuneShort and CurrentRecoveryTime as well as
        % AgentCharacteristics.ImmunityShort and
        % AgentCharacteristics.RecoveryTime: 
        %%%%%%CurrentImmuneShort = AgentCharacteristics.ImmunityShort;
        %%%%%%CurrentRecoveryTime = AgentCharacteristics.RecoveryTime;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% RECOVERY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

        % Find index of agents that are infected   
        IAgentsInf = find(CurrentInfections > 0);

        %if isempty(IAgentsInf)
            %fprintf('Pathogen extinct')
        %    break
        %end

        RecoveryProbability = parameters.Precovery * ones(size(CurrentInfections));
 
        % Account for short ss-immunity in RecoveryProbability:
        if parameters.SigmaShort > 0
            RecoveryProbability(CurrentImmuneShort > 0) = parameters.PrecoveryImmuneSSS;
        end

        % Account for long ss-immunity in RecoveryProbability:
        if parameters.SigmaLong < 1
            if parameters.SigmaLong > 0
                RecoveryProbability(CurrentImmuneLong > 0) = parameters.PrecoveryImmuneSSL;
            end
        end

        % Account for short cs-immunity in RecoveryProbability:
        if parameters.OmegaShort > 0
            CurrentImmunityToAnyStrain = repmat(any(CurrentImmuneShort > 0,2),1,parameters.MaxStrains);
            RecoveryProbability(CurrentImmuneShort == 0 & CurrentImmuneLong == 0 & CurrentImmunityToAnyStrain == 1) = ...
                parameters.PrecoveryImmuneCSS; 
        end

        % Account for long cs-immunity in RecoveryProbability:
        if parameters.OmegaLong > 0
            CurrentImmunityToAnyStrain = repmat(any(CurrentImmuneLong > 0,2),1,parameters.MaxStrains);
            RecoveryProbability(CurrentImmuneLong == 0 & CurrentImmuneShort == 0 & CurrentImmunityToAnyStrain == 1) = ...
                parameters.PrecoveryImmuneCSL;
        end
      
        % Remove index of agents from IAgentsInf that do not recover    
        IAgentsInf(rand(size(IAgentsInf)) < (1 - RecoveryProbability(IAgentsInf))) = [];

        % Update infection and immune status
        AgentCharacteristics.Infections(IAgentsInf) = 0;
        AgentCharacteristics.ImmunityShort(IAgentsInf) = 1 + ...
            CurrentImmuneShort(IAgentsInf);
        AgentCharacteristics.RecoveryTime(IAgentsInf) = CurrentTime;

        % Check if agents have obtained permanent immunity
        IAgentsImm = find(AgentCharacteristics.ImmunityShort == parameters.NumberRepeatInfections);

        % Update agents that have obtained permanent immunity  
        AgentCharacteristics.ImmunityLong(IAgentsImm) = 1;
        AgentCharacteristics.ImmunityShort(IAgentsImm) = 0;
        AgentCharacteristics.RecoveryTime(IAgentsImm) = 0;  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TRANSMISSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        % Find row indices of agents that are infected. And 
        % store in RowIndexInfected.
        % Start off with all agents and remove those which are not infected
        RowIndexInfected = KeepIndexAllAgents'; % all agents
        G = sum(CurrentInfections,2);
        RowIndexInfected(G==0)=[];

        if ~isempty(RowIndexInfected)  

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Determine susceptibility of each agent to each strain:    
            % InfectionProbability stores the probability that each agent
            % will contract each strain following a contact

            % Infection status of agent:  CurrentInfections
            % Number of co-infections of each agent:
            TotalInfections = sum(CurrentInfections,2);

            % Immunity status of agent: CurrentImmune
            % Account for resistance to co-infection in Ptransmission1
            Ptransmission1 = parameters.Ptransmission .* ...
                                (1 - TotalInfections / parameters.CCC) .^ (parameters.coinfection);  
            Ptransmission1 = repmat(Ptransmission1,1,parameters.MaxStrains);
            InfectionProbability = Ptransmission1;

            % When ss-immunity is 100% effective, and immunity 
            % reduces duration of infection, this is equivalent 
            % to being 100% protected against infection: 

            if parameters.SigmaLong == 1
                InfectionProbability(CurrentImmuneLong > 0) = ...
                 Ptransmission1(CurrentImmuneLong > 0) * 0;
            end

             % For each agent that is infected
             for j = 1 : length(RowIndexInfected)

                % Here are the strains infecting the agent:
                InfectingStrainsAgent = find(CurrentInfections(RowIndexInfected(j),:));

                % Determine number of contacts with other agents
                X = ContactRand(countCR,1);
                countCR = countCR + 1;
                if countCR > length(ContactRand)
                    ContactRand = poissrnd(parameters.Ncontacts,1e6,1);
                    countCR = 1;
                end

                % If the agent makes contact with other agents 
                if X > 0
                    if parameters.PopSize > 1
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Find out which agents are contacted. Here,
                        % IndexOfContacts = orignal row indices of contacted agents
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Sample X agents with replacement (using 
                        % pregenerated random numbers)from all agents in the
                        % population
                        IndexOfContacts = ceil(parameters.PopSize * ...
                            SamplingContactsRand(countSCR:countSCR+X-1,1));
                        % Update counter for pregenerated random numbers:
                        countSCR = countSCR+X;
                        % Regenerate random numbers if they are all used up:
                        if countSCR > length(SamplingContactsRand) - max(parameters.PopSize)
                            SamplingContactsRand = rand(1e6,1);
                            countSCR = 1;
                        end
                        % if this sample includes the infected agent of interest,
                        % remove these and resample:
                        if ~isempty(IndexOfContacts(IndexOfContacts==RowIndexInfected(j)))
                            sampleextra = ...
                                datasample(KeepIndexAllAgents(KeepIndexAllAgents~=RowIndexInfected(j)),...
                                        length(IndexOfContacts(IndexOfContacts==RowIndexInfected(j))));

                            %a = sampleextra(:)
                            %b = IndexOfContacts(IndexOfContacts~=RowIndexInfected(j))
                            IndexOfContacts = [IndexOfContacts(IndexOfContacts~=RowIndexInfected(j));...
                                                sampleextra(:)];
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

                        InfectingStrains = ceil(length(InfectingStrainsAgent)*SamplingContactsRand(countSCR:countSCR+X-1,1));

                        if length(InfectingStrainsAgent) == 1
                            InfectingStrains = InfectingStrainsAgent(InfectingStrains);
                        else
                            InfectingStrains = InfectingStrainsAgent(InfectingStrains)';
                        end

                        % Update counter for pregenerated random numbers:
                        countSCR = countSCR+X;

                        % Regenerate random numbers if they are all used up:
                        if countSCR > length(SamplingContactsRand) - parameters.PopSize
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
                        if countSCR > length(SamplingContactsRand) - parameters.PopSize
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
                            %idx = sub2ind(size(AgentCharacteristics.Infections), IndexOfContacts, InfectingStrains); <- SLOW
                            idx = IndexOfContacts + (InfectingStrains - 1) * size(AgentCharacteristics.Infections,1); % <- FAST
                            AgentCharacteristics.Infections(idx) = CurrentInfections(idx) + 1;

                            for l = 1:length(idx)
                                if CurrentImmuneShort(idx(l))==1
                                    AgentCharacteristics.RecoveryTime(idx(l)) = NaN;
                                end                          
                            end

                            % Make infected susceptibles not susceptible for next 
                            % infected agent of interest:
                            InfectionProbability(IndexOfContacts,:) = ...
                                zeros(length(IndexOfContacts),parameters.MaxStrains);

                        end

                    end
                end

            end

        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% AGING, BIRTH & DEATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for k = 1 : parameters.NPopulations   
            % Aging:
            AgentCharacteristics.Age = parameters.dt_years + ...
                            AgentCharacteristics.Age;

            % Death and Birth: remove agents that are older than 
            % AgeDeath, replace with susceptible newborn (age 0.001)
            D = find(AgentCharacteristics.Age > parameters.AgeDeath);
            %store_deaths = [store_deaths; D];
            AgentCharacteristics.Infections(D,:) = 0;
            AgentCharacteristics.ImmunityShort(D,:) = 0;
            AgentCharacteristics.ImmunityLong(D,:) = 0;
            AgentCharacteristics.VaccinationStatus(D,1) = 0;
            AgentCharacteristics.RecoveryTime(D,:) = 0;
            AgentCharacteristics.Age(D) = 0.001;
            
            % Routine vaccination
            if intervention == 2

                if i > parameters.iIntervention

                    D = find(AgentCharacteristics.Age > parameters.AgeOfVacccination & ...
                        AgentCharacteristics.Age <= parameters.AgeOfVacccination + parameters.dt_years);
                    Dt = randperm(length(D),round(parameters.RoutineVaccinationCoverage * length(D)));
                    D = D(Dt);
                    AgentCharacteristics.VaccinationStatus(D,1) = ones(length(D),1);
                    AgentCharacteristics.ImmunityLong(D,parameters.VaccineStrains) = ...
                        ones(length(D),parameters.NumberOfVaccineStrains);
                    AgentCharacteristics.ImmunityShort(D,parameters.VaccineStrains) = ...
                        zeros(length(D),parameters.NumberOfVaccineStrains);
                    AgentCharacteristics.Infections(D,parameters.VaccineStrains) = ...
                        zeros(length(D),parameters.NumberOfVaccineStrains);

                end
            end
        end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% IMMIGRATION FROM OUTSIDE ALL POPULATIONS %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Remove agents randomly, replace with agents with
        % random age, new (or extinct) strain infection, and no immune history
        
        % Determine number of immigrants from pregenerated random numbers:
        NumMig = MRRand(countMR);
        % Update counter for pregenerated random numbers:
        countMR = countMR + 1;
        % Regenerate random numbers if they are all used up:
        if countMR > 1e6
            MRRand = poissrnd(parameters.ImmigrationRatePerCapita * ...
                            parameters.PopSize,1e6,1);
            countMR = 1;
        end
        
        % Determine infection status of immigrants:
        Infected_migrants = (rand(NumMig,1) < parameters.prevalence_in_migrants);

        % Number of infected immigrants:
        NumIMig = sum(Infected_migrants);
        
        % Determine strains infecting immigrants:
%         TotalNSPresent = SummaryStatistics.NumberVaccineStrainsPresent(1,i) + ...
%                             SummaryStatistics.NumberNonVaccineStrainsPresent(1,i);
%         ProbStrainIsVS = (SummaryStatistics.NumberVaccineStrainsPresent(1,i) + 1) / ...
%                             (TotalNSPresent + 2);
%         ProbStrainIsVS(ProbStrainIsVS > 1) = 1;
%         NumberMigrantStrainsVS = binornd(NumIMig,ProbStrainIsVS);
%         NumberMigrantStrainsNVS = NumIMig - NumberMigrantStrainsVS;
%         SelectVS = randperm(parameters.NumberOfVaccineStrains,NumberMigrantStrainsVS);
%         SelectNVS = randperm(parameters.MaxStrains - parameters.NumberOfVaccineStrains, ...
%                         NumberMigrantStrainsNVS);
%         SelectedVS = parameters.VaccineStrains(SelectVS);
%         SelectedNVS = parameters.NVaccineStrains(SelectNVS);  
%         
%         MigrantStrains = [SelectedVS(:)' SelectedNVS(:)'];
        
        SelectedStrains = randperm(parameters.MaxStrains, NumIMig);
        MigrantStrains = SelectedStrains(:)';
        
        countM = 1; % counts number of immigrants into all populations
        countIM = 1; % counts number of infected immigrants into all populations
            
        % Agents leaving population:
        D = randperm(parameters.PopSize,NumMig)';
        %store_migrants = [store_migrants; D];

        % Update infection, immune, age status of the population:
        for h = 1 : NumMig

            % immigrants have up to one infection:
            AgentCharacteristics.Infections(D(h),:) = zeros(1,parameters.MaxStrains);
            if Infected_migrants(countM)>0
                AgentCharacteristics.Infections(D(h),MigrantStrains(countIM)) = 1;                  
                countIM = countIM + 1;                    
            end

            AgentCharacteristics.ImmunityLong(D(h),:) = zeros(1,parameters.MaxStrains);

            % Immigrants have same vaccination status and long immune status
            % as that of an agent randomly selected from the population
            randagent = ceil(rand * parameters.PopSize);
            AgentCharacteristics.ImmunityLong(D(h),:) = CurrentImmuneLong(randagent,:);
            AgentCharacteristics.VaccinationStatus(D(h),1) = CurrentVaccineCoverage(randagent,1);

            % immigrants have no short term immunity:
            AgentCharacteristics.ImmunityShort(D(h),:) = zeros(1,parameters.MaxStrains);          
            AgentCharacteristics.RecoveryTime(D(h),:) = zeros(1,parameters.MaxStrains);

            % immigrants have random age:
            AgentCharacteristics.Age(D(h),end) = rand * parameters.AgeDeath;
            countM = countM + 1;

        end
                   

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PREVALENCE and COINFECTION DISTRIBUTIONS %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

        SummaryStatistics.SingleStrainPrevTime(:,i + 1) = sum(AgentCharacteristics.Infections,1)';
        SummaryStatistics.SingleStrainLongImmunityTime(:,i + 1) = sum(AgentCharacteristics.ImmunityLong,1)';
        
        [counts,K] = calculate_coinfection_distribution(AgentCharacteristics.Infections); 

        SummaryStatistics.CoinfectionDistributionTime(K,i + 1) = counts;

        [counts,K] = calculate_coinfection_distribution(AgentCharacteristics.ImmunityLong);  

        SummaryStatistics.LongImmunityDistributionTime(K,i + 1) = counts;

        SummaryStatistics.AvgLongImmunity(1,i+1) = ...
            calculate_avg_immunity(SummaryStatistics.LongImmunityDistributionTime(:,i + 1)) ./ parameters.PopSize;

        SummaryStatistics.PopulationSizeTime(1,i + 1) = parameters.PopSize;

        SummaryStatistics.AvgVaccinationCoverage(1,i + 1) =  mean(AgentCharacteristics.VaccinationStatus);

        [SummaryStatistics.PrevalenceVaccinated(1,i + 1),...
        SummaryStatistics.PrevalenceUnvaccinated(1,i + 1),...
        SummaryStatistics.NumberVaccinated(1,i + 1),...
        SummaryStatistics.NumberUnvaccinated(1,i + 1)] = ...
            calculate_prev_vacc(AgentCharacteristics.Infections,...
                                AgentCharacteristics.VaccinationStatus);

        [counts,K] = calculate_coinfection_distribution(AgentCharacteristics.ImmunityLong(:,parameters.VaccineStrains));  
        SummaryStatistics.LongImmunityDistributionTimeVS(K,i + 1) = counts;

        [counts,K] = calculate_coinfection_distribution(AgentCharacteristics.ImmunityLong(:,parameters.NVaccineStrains));  
        SummaryStatistics.LongImmunityDistributionTimeNVS(K,i + 1) = counts;

        SummaryStatistics.AvgLongImmunityVS(1,i + 1) =  calculate_avg_immunity(SummaryStatistics.LongImmunityDistributionTimeVS(:,i + 1)) ./ parameters.PopSize;
        SummaryStatistics.AvgLongImmunityNVS(1,i + 1) =  calculate_avg_immunity(SummaryStatistics.LongImmunityDistributionTimeNVS(:,i + 1)) ./ parameters.PopSize;

        [SummaryStatistics.NumberVaccineStrainsPresent(1,i + 1),...
        SummaryStatistics.NumberNonVaccineStrainsPresent(1,i + 1)] = ...
            calculate_number_strain_present(AgentCharacteristics.Infections,...
            parameters.VaccineStrains, parameters.NVaccineStrains);

        SummaryStatistics.AgeDistributionPrevalence(:,i + 1) = ...
            calculate_age_dist_prev(AgentCharacteristics.Infections,AgentCharacteristics.Age,parameters.AgeCategories);   

    end
    
    SSPrev = SummaryStatistics.SingleStrainPrevTime;
    AgentsInfectedByKStrains = SummaryStatistics.CoinfectionDistributionTime;

    function [counts, K] = calculate_coinfection_distribution(Infections)

        if sum(sum(Infections))>1
            counts = sum(Infections,2);
            counts(counts==0)=[];
            [counts,K] = hist(counts,unique(counts)); % Counts number of agents (counts) with K infections
            if length(counts(counts>0))==1
                counts = counts(counts>0);
                K = K(counts>0);
            end
        elseif sum(sum(Infections))==1
            counts = 1; 
            K = 1;
        else
            counts = 0; 
            K = 1;
        end

    end

    function avg = calculate_avg_immunity(Immunity)
        ms = length(Immunity);
        ms = (1:1:ms)';
        avg = sum(ms .* Immunity,1);
    end

    function [VacPrev, UnVacPrev, NumVac, NumUVac] = calculate_prev_vacc(ACI,ACVS)
        tempAC = sum(ACI,2);
        NumVac = sum(ACVS);
        NumUVac = size(ACI,1) - NumVac;
        NumInfVac = length(tempAC(ACVS==1 & tempAC>0));
        NumInfUVac = length(tempAC(ACVS==0 & tempAC>0));
        if NumVac > 0
            VacPrev = NumInfVac / NumVac * 100; 
        else
            VacPrev = 0;
        end
        if NumUVac > 0
            UnVacPrev = NumInfUVac / NumUVac * 100; 
        else
            UnVacPrev = 0;
        end
    end

    function [NumVacStrainP, NumNonVacStrainP] = calculate_number_strain_present(ACI,VS,NVS)
        
        ACI = sum(ACI,1);
        ACI(ACI>0) = 1;
        NumVacStrainP = sum(ACI(1,VS),2);
        NumNonVacStrainP = sum(ACI(1,NVS),2);
        
    end

    function AgeDistPrev = calculate_age_dist_prev(ACI,Age,AgeCat)
        ACI = sum(ACI,2);
        ACI(ACI>0) = 1;
        agePrev = Age.*ACI;
        agePrev(agePrev==0) = [];
        AgeDist = histcounts(Age,AgeCat);
        AgeDistPrev = histcounts(agePrev,AgeCat);
        AgeDistPrev = AgeDistPrev(:) ./ AgeDist(:);       
    end

end

