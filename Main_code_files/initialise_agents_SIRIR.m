function AgentCharacteristics = initialise_agents_SIRIR(parameters)
    
    Infections = zeros(parameters.PopSize, parameters.MaxStrains);
    Age = rand(parameters.PopSize,1) * parameters.AgeDeath;
    VS = zeros(parameters.PopSize,1);

    % AgentCharacteristics.XXX(i,j): is the XXX status of agent i with 
    % respect to the jth strain. Specifically,  
    % AgentCharacteristics.Infections(i,j): 0 = susceptible, 
    % n = infected by n copies of strain j.  
    % AgentCharacteristics.RecoveryTime(i,j): time that host i obtained 
    % short-term immunity to strain j (if 0, currently does not have 
    % short-term immunity).
    % AgentCharacteristics.ImmunityShort(i,j): 0 = no short term immunity, 
    % >1 = short term immunity to strain j.  Once this reaches 
    % parameters.NumberRepeatInfections, then long-term immunity is
    % obtained.
    % AgentCharacteristics.ImmunityLong(i,j): 0 = no life-long immunity, 
    % 1 = life-long immunity to strain j.    
    % AgentCharacteristics.Age(i): is the age of agent i 
    % and is sampled from a uniform distribtion
    
    AgentCharacteristics.Infections = Infections;
    AgentCharacteristics.RecoveryTime = Infections;
    AgentCharacteristics.ImmunityShort = Infections;
    AgentCharacteristics.ImmunityLong = Infections;
    AgentCharacteristics.Age = Age;
    AgentCharacteristics.VaccinationStatus = VS;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIAL INFECTION STATUS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Randomly select NI0perstrain agents to be infected by each strain in 
    % population, and randomly select RI0perstrain to have life-long 
    % immunity to each strain in population.
        
    % Assign an index to each agent
    AgentIndexInf = 1 : 1 : parameters.PopSize; 
    AgentIndexImm = AgentIndexInf;

    % Determine initial strains in population (immunity can be to any 
    % strain in all populations)
    Initialstrains = randperm(parameters.MaxStrains,parameters.Nstrains);

    % For each strain 
    if ~isempty(Initialstrains)
        for i = 1 : parameters.Nstrains

            strain = Initialstrains(i);

            % Choose agents to be infected
            InfectedAgents = randperm(length(AgentIndexInf),parameters.NI0perstrain); 

            % Update agent's characteristics
            AgentCharacteristics.Infections(AgentIndexInf(InfectedAgents),strain) = ...
                AgentCharacteristics.Infections(AgentIndexInf(InfectedAgents),strain) + 1;         

            % Sample without replacement 
            %AgentIndexInf(InfectedAgents) = [];

        end
    end

    % For each strain 
    for i = 1 : parameters.MaxStrains 

        % Choose agents to be immune
        ImmuneAgents = randperm(length(AgentIndexImm),parameters.NR0perstrain); 

        % Update agent's characteristics         
        AgentCharacteristics.ImmunityLong(AgentIndexImm(ImmuneAgents),i) = 1; 

        % Sample without replacement 
        % AgentIndexImm(ImmuneAgents) = [];

    end


end