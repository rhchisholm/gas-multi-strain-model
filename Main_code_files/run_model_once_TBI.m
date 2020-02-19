function [TimeBetweenInfectionSS,AvgPrev,AvgDiv,PeriodDiv] = run_model_once_TBI(R0,DI)

    sigma = 1;
    omega = 0.1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Duration of simulation (years)
    DurationSimulation = 100;
    % Number of strains
    Nstrains = 42;
    % Total number of hosts
    Nagents = 2500; 
    % Age that all hosts die (years)
    AgeDeath = 71; 
    % Basic Reproduction Number 
    BasicReproductionNumber = R0 ;  
    % Duration of immunity (weeks)
    Dimmunity = DI * 52.14; 
    % Resistance to co-infection
    x = 10;
    % Migration rate per week per population
    alpha = 3;  % daily per capita rate is alpha / 7 / Nagents
    % Number of contacts per week
    Cperweek = 34.53;

    params = double([DurationSimulation, Nstrains, Dimmunity, ...      
                sigma, omega, x, ...              
                Cperweek, Nagents, alpha, ... 
                AgeDeath, BasicReproductionNumber]);                               

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialise agents %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [AgentCharacteristics, ImmuneStatus, ~] = initialise_agents(params);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Simulate transmission %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [SSPrev,AgentsInfectedByKStrains,TimeBetweenInfectionSS] = ...
                    simulator_TBI(AgentCharacteristics, ImmuneStatus, params, 0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate summary stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time = linspace(0,DurationSimulation,size(SSPrev,2));

    Diversity = div(SSPrev);
    AvgDiv = Diversity(time>30);
    Prevalence = sum(AgentsInfectedByKStrains,1)*100/Nagents;
    AvgPrev = Prevalence(time>30);
    %Diversity rolling average over 10 years from prevalence data
    for i = 1:length(time)
        if i <= 365*10
            timesteps = 1:1:i;
        else
            timesteps = i-365*10+1:1:i;
        end
        PrevD(:,i)=sum(SSPrev(:,timesteps),2);
    end
    DivAvg2 = div(PrevD);
    PeriodDiv = DivAvg2(time>30);
    
    % Calculates diversity given number of infections of each strain
    function D = div(SSP)

        SSP1 = SSP - 1;
        N = sum(SSP);
        D = N .* (N - 1);
        sumSSP = sum(SSP .* SSP1);
        D = D ./ sumSSP;
        D(D == Inf) = N(D == Inf); 
        D(isnan(D)) = 0;

    end

end

