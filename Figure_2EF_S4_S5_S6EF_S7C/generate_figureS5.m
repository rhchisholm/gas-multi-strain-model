% Generates Figure S5 

clear all
close all

set(0,'DefaultTextFontName','Arial')
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultAxesFontName','Arial')

if ~exist('figure_files', 'dir')
    mkdir('figure_files')
end

% Varying R0 figure

% w = 1: x = 10
% w = 2: x = 100

for w = 1 : 2
    
    %%%%%%%%%%%%%
    % LOAD DATA %
    %%%%%%%%%%%%%
    
    if w == 1 
        load('mat_files/figure2E.mat');
        load('mat_files/figureS5A.mat');
    elseif w == 2
        load('mat_files/figureS4A.mat');
        load('mat_files/figureS5C.mat');
    end
    
    R0 = 0.0301 / (1/71/52.14 + 1/2 + params(9) / params(8));
    
    TotalAgentsInfected = zeros(length(Dimmunity),length(Contacts),NumberSimulations);
    TotalPrev = zeros(length(Dimmunity),length(Contacts),NumberSimulations);
    Diversity = zeros(length(Dimmunity),length(Contacts),NumberSimulations);
    
    for i = 1 : length(Dimmunity)
            
        for j = 1 : length(Contacts)

            SSPrevTimePoint(:,:,i,j) = squeeze(TimePrev(i,j,:,:));
            NumAgentsInfected(:,:,i,j) = TimeAgentsInfectedByKStrains(i,j,:,:);

            for k = 1 : NumberSimulations

                TotalPrev(i,j,k) = sum(SSPrevTimePoint(:,k,i,j),1);
                TotalAgentsInfected(i,j,k) = sum(NumAgentsInfected(:,k,i,j),1);
                Diversity(i,j,k) = div(squeeze(SSPrevTimePoint(:,k,i,j)));
            
            end

        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE DIVERSITY, PREVALENCE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
    % Calculate mean final diversity
    z = mean(Diversity,3);
    z = squeeze(z);

    % Calculate mean final prevalence of each simulation
    AvgTotalAgents = mean(TotalAgentsInfected,3);
    Prevalence = AvgTotalAgents / Nagents;

    figure(w)

    for i = 1 : length(Dimmunity)
        if i == 2 || i == 3
            hl = scatter(z(i,:),Prevalence(i,:)*100,200,Contacts * R0,'filled');
        else
            hl = scatter(z(i,:),Prevalence(i,:)*100,200,Contacts * R0);
        end
        Markers = {'o','s','p','d'};
        set(hl, 'marker', strcat(Markers{i}))
        set(gca, 'YScale', 'log')
        hold on
    end
    
    ylabel('Mean endemic prevalence (%)')
    xlabel('Mean endemic diversity')
    axis([0 NumberInitialStrains 0.2 100])
    h = colorbar;
    ylabel(h, {'Basic reproduction number (R_0)'})
    caxis([1 10.5])  
    legend('w = 3 weeks','w = 19 weeks','w=104 weeks','Location','southeast')
    
    if w == 1
        savefig('figure_files/figureS5A.fig')
        saveas(gcf,'figure_files/figureS5A','epsc')
    elseif w == 2
        savefig('figure_files/figureS5C.fig')
        saveas(gcf,'figure_files/figureS5C','epsc')
    end   
end


% Varying N figure

% w = 1: x = 10
% w = 2: x = 100

N = [200 1000 2000 3500 5000 10000];

for w = 1 : 2

    for di = 1:3
    %%%%%%%%%%%%%
    % LOAD DATA %
    %%%%%%%%%%%%%
        if w == 1
            load('mat_files/figureS5B.mat');
        elseif w == 2            
            load('mat_files/figureS5D.mat');
        end
    
        TotalAgentsInfected = zeros(length(N),NumberSimulations);
        TotalPrev = zeros(length(N),NumberSimulations);
        Diversity = zeros(length(N),NumberSimulations);

        for j = 1 : length(N)

            SSPrevTimePoint(:,:,j) = squeeze(TimePrev(1,j,:,:));
            NumAgentsInfected(:,:,j) = TimeAgentsInfectedByKStrains(1,j,:,:);

            for k = 1 : NumberSimulations

                TotalPrev(j,k) = sum(SSPrevTimePoint(:,k,j),1);
                TotalAgentsInfected(j,k) = sum(NumAgentsInfected(:,k,j),1);
                Diversity(j,k) = div(squeeze(SSPrevTimePoint(:,k,j)));

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CALCULATE DIVERSITY, PREVALENCE %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Calculate mean final diversity
        z = mean(Diversity,2);
        z = squeeze(z);

        % Calculate mean final prevalence of each simulation
        AvgTotalAgents = mean(TotalAgentsInfected,2);
        Prevalence = AvgTotalAgents ./ N';

        figure(3 + w)

        if di == 2 || di == 3
            hl = scatter(z,Prevalence*100,200,N,'filled');
        else
            hl = scatter(z,Prevalence*100,200,N);
        end

        Markers = {'o','s','p','d'};
        set(hl, 'marker', strcat(Markers{di}))
        set(gca, 'YScale', 'log')
        hold on
    end

    ylabel('Mean endemic prevalence (%)')
    xlabel('Mean endemic diversity')
    axis([0 NumberInitialStrains 0.2 100])
    h = colorbar;
    ylabel(h, {'Population size (N)'})
    caxis([100 10000])
    set(gca,'ColorScale','log')
    legend('w = 3 weeks','w = 19 weeks','w=104 weeks','Location','southeast')
    
    if w == 1
        savefig('figure_files/figureS5B.fig')
        saveas(gcf,'figure_files/figureS5B','epsc')
    elseif w ==2
        savefig('figure_files/figureS5D.fig')
        saveas(gcf,'figure_files/figureS5D','epsc')
    end
end


function D = div(SSP)

    SSP1 = SSP - 1;
    N = sum(SSP);
    D = N .* (N - 1);
    sumSSP = sum(SSP .* SSP1);
    D = D ./ sumSSP;
    D(D == Inf) = N(D == Inf);
    D(isnan(D)) = 0;

end
