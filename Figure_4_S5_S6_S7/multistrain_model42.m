function SumStats = multistrain_model42(R0,DI,sigma,omega,x,NumMigrationEvents,burnin_time)
    
    rng(1)
    
    %burnin_time = 20.0;
    
    if DI < 0
    	DI = 0;
    end

    params = double([burnin_time, 42, DI*52.14, ...      % Nstrains
                sigma, omega, x, ...              
                30, 2500, NumMigrationEvents, ... % Cperweek, Na
                71, R0]);                          % AgeDeath
            
    EnrolledPeople = 548;

    Consultations = [27,21,42,51,36,69,...
                    122,149,172,170,142,147,...
                    40,193,183,211,190,182,...
                    199,130,191,188,161]; % number of consultations at each time point

    time_obs = [1,32,62,93,123,154,...
                185,214,245,275,306,336,...
                367,398,428,459,489,520,...
                551,579,610,640,671]; % days of consultations
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Burn in to reach enedmic equilibrium %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %NumMigrationEvents
    
    [AgentCharacteristics, ImmuneStatus, ~] = initialise_agents(params);

    [~, ~, ~, ~, AgentCharacteristics, ImmuneStatus, params] = ...
        simulator2(AgentCharacteristics, ImmuneStatus, params, 1, Consultations, time_obs);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate synthetic data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %size(ImmuneStatus)
    %size(AgentCharacteristics)
    
    params(1) = 2;

    [~, ~, Observed_AC_time, ~, ~, ~, params] = ...
            simulator2(AgentCharacteristics, ImmuneStatus, params, 0, Consultations, time_obs);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Sample synthetic data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This vector stores the index of enrolled agents 
    Index_enrolled = randperm(params(8),EnrolledPeople);

    % This will store the infection status of observed agents at each
    % time point
    AC_obs = zeros(size(Observed_AC_time));

    % At each time point:
    for j = 1 : length(time_obs)

        % This is the number of consulations:
        C = Consultations(j);
        % Randomly choose consultations from enrolled population:
        Index_agents_observed = datasample(Index_enrolled,C,'Replace',false);
        % And store their infection status at the jth time point
        AC_obs(1:C,:,j) = Observed_AC_time(Index_agents_observed,:,j);

    end

    % Store the total number of each strain observed at each time point
    % in the ith sampling process
    SSPrev_obs = squeeze(sum(AC_obs,1));
    
    %size(SSPrev_obs)

    % Calculate total prevalence and diversity
    for j = 1 : length(time_obs)

        [IAgent, ~] = find(squeeze(AC_obs(:,:,j)) > 0);
        Prevalence_obs(j) = length(unique(IAgent)) / ...
                                Consultations(j) * 100;
        Diversity_obs(j) = div(squeeze(SSPrev_obs(:,j)));

    end
    
    PD = [Prevalence_obs; Diversity_obs];
    
    % Summary statistics
    
    AvgTimeObsStr = SSPrev_obs;
    AvgTimeObsStr(SSPrev_obs>0)=1;
    AvgTimeObsStr = sum(AvgTimeObsStr,2);
    MaxTimeObsStr = max(AvgTimeObsStr);
    AvgTimeObsStr = mean(AvgTimeObsStr(AvgTimeObsStr>0));
    NumStrainsObs = sum(SSPrev_obs,2);
    NumStrainsObs = length(NumStrainsObs(NumStrainsObs>0));
    AvgTimeRepeatInf = mean(timerepeat(SSPrev_obs));
    VarTimeRepeatInf = var(timerepeat(SSPrev_obs));
    

    AvgPrev = mean(Prevalence_obs);
    AvgDiv = mean(Diversity_obs);
    MaxAbundance = max(max(SSPrev_obs));
    VarPrev = var(Prevalence_obs);
    VarDiv = var(Diversity_obs);
    AvgNPMI = MI(SSPrev_obs);
    
    DivAllIsolates = div(sum(SSPrev_obs,2));
    

    SumStats = [AvgPrev;AvgDiv;MaxAbundance;AvgTimeObsStr;...
                MaxTimeObsStr;NumStrainsObs;VarPrev;VarDiv;...
                AvgTimeRepeatInf; VarTimeRepeatInf; AvgNPMI;...
                DivAllIsolates];
    
    SumStats(isnan(SumStats)) = 0;

    
    % Calculates diversity given number of infections of each strain
    function D = div(SSP)

        SSP1 = SSP - 1;
        N = sum(SSP);
        D = N .* (N - 1);
        sumSSP = sum(SSP .* SSP1);
        D = D ./ sumSSP;
        D(D == Inf) = N(D == Inf); % fixed calculation when there are one of each strain (n/0)
        D(isnan(D)) = 0;
        %D(D > 0) = D(D > 0) ./ sumSSP(D > 0);
        %D(D == Inf) = 0;

    end

    function E = timerepeat(SSPrev)
        A = find(SSPrev' > 0);
        E = [];
        for i = 1:size(SSPrev,1)
            vec = A(A <= i * size(SSPrev,2) & A >= ((i-1)*size(SSPrev,2) + 1));
            vec = diff(vec);
            vec = vec(vec > 1);
            vec = vec - 1;
            E = [E; vec];
        end     
    end
    
    % Calculates normalised pointwise mutual information
    function F = MI(SSP)
    	TotalO = sum(sum(SSP));
    	PX = sum(SSP,2)/TotalO;
    	PXY = zeros(size(SSP,1));
    	NPMI = PXY;
    	for i = 1:size(SSP,1)
    		for jj = 1:size(SSP,1)
    			if jj < i
    				PXY(i,jj) = sum(min(SSP(i,:),SSP(jj,:)))/TotalO;
    			end
    		end
    	end
    
    	for i = 1:size(SSP,1)
    		for jj = 1:size(SSP,1)
    			if jj < i
    				if (PX(i)+PX(i)) > 0
    					if PXY(i,jj) > 0
    						NPMI(i,jj) = (-log(PX(i))-log(PX(jj))+log(PXY(i,jj)))/(-log(PXY(i,jj)));
    					else
    						NPMI(i,jj) = -1;
    					end 
    				else
    					NPMI(i,jj) = 0;
    				end
    			end
    		end
    	end
    
    	F = sum(sum(NPMI))/((size(SSP,1)^2-size(SSP,1))/2);

    end


end
