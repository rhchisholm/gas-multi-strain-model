function [AgentCharacteristics,ImmuneStatus,time] = initialise_agents(params)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialise agents %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %parameters;

    [Nagents , Nstrains, Nst, AgeDeath, NI0perstrain, NR0perstrain,...
        ~, ~, ~, ~, ...
        ~, ~, ~, ~, ...
        ~, ~, ~, ...
        time, ~, ~] = parameters(params);

    % AgentCharacteristics(i,j): jth characteristic of agent i
    % forall j in [1, Nstrains], AgentCharacteristics(i,j) 
    % reflects the status of agent i with respect to the jth  
    % strain: 0 = susceptible, n = infected by n copies of
    % strain j.
    AgentCharacteristics = zeros(Nagents , Nstrains + 1); 

    % AgentCharacteristics(i,end) is the age of agent i. Sample 
    % ages from uniform distribtion:
    Ages = rand(1,Nagents) * AgeDeath;
    AgentCharacteristics(:,end) = Ages';
    clear Ages

    % ImmuneStatus(i,j): immune status of agent i wrt 
    % strain j, forall j in [1, Nstrains]: 
    % 0 = no strain-specific immunity
    % 1 = strain-specific immunity of strength StrengthImmunity
    ImmuneStatus = zeros(Nagents , Nstrains); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIAL INFECTION STATUS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Randomly select NI0perstrain agents to be infected by each 
    % strain.

    % Assign an index to each agent
    AgentIndexInf = 1 : 1 : Nagents; 
    AgentIndexImm = AgentIndexInf;
    
    if Nagents<NI0perstrain*Nst
        NI0perstrain = 4;
        NR0perstrain = 4;
    end
    
    % For each strain 
    for i = 1 : Nst 

        % Choose agents to be infected
        InfectedAgents = randperm(length(AgentIndexInf),NI0perstrain); 

        % Choose agents to be immune
        ImmuneAgents = randperm(length(AgentIndexImm),NR0perstrain); 

        % Update agent's characteristics
        AgentCharacteristics(AgentIndexInf(InfectedAgents),i) = 1;        
        ImmuneStatus(AgentIndexImm(ImmuneAgents),i) = 1; 

        % Sample without replacement 
        AgentIndexInf(InfectedAgents) = [];
        AgentIndexImm(ImmuneAgents) = [];

    end




    clear AgentIndexInf AgentIndexImm InfectedAgents ImmuneAgents

end