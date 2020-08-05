% Generates Figure 2E-F, S4A-B and S6E-F 

clear all
close all

set(0,'DefaultTextFontName','Arial')
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultAxesFontName','Arial')

if ~exist('figure_files', 'dir')
    mkdir('figure_files')
end

% Varying R0 figures: 2E, S4A, S6E

% w = 1: x = 10, extra_ci_effect = 0;
% w = 2: x = 100, extra_ci_effect = 0;
% w = 3: x = 10, extra_ci_effect = 1;

for w = 1 : 3
    
    %%%%%%%%%%%%%
    % LOAD DATA %
    %%%%%%%%%%%%%
    
    if w == 1 
        load('mat_files/figure2E.mat');
    elseif w == 2
        load('mat_files/figureS4A.mat');
    else
        load('mat_files/figureS6E.mat');
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
    legend('1/\theta = 0 years','1/\theta = 6 months','1/\theta = 5 years','1/\theta = lifelong','Location','southeast')
    
    if w == 1
        savefig('figure_files/figure2E.fig')
        saveas(gcf,'figure_files/figure2E','epsc')
    elseif w == 2
        savefig('figure_files/figureS4A.fig')
        saveas(gcf,'figure_files/figureS4A','epsc')
    else
        savefig('figure_files/figureS6E.fig')
        saveas(gcf,'figure_files/figureS6E','epsc')

    end   
end


% Varying N figures: 2F, S4B, S6F

% w = 1: x = 10, extra_ci_effect = 0;
% w = 2: x = 100, extra_ci_effect = 0;
% w = 3: x = 10, extra_ci_effect = 1;

N = [200 1000 2000 3500 5000 10000];

for w = 1 : 3

    for di = 1:4
    %%%%%%%%%%%%%
    % LOAD DATA %
    %%%%%%%%%%%%%
        if w == 1
            if di == 1
                load('mat_files/figure2F_1.mat');
            elseif di ==2
                load('mat_files/figure2F_2.mat');
            elseif di == 3
                load('mat_files/figure2F_3.mat');
            else
                load('mat_files/figure2F_4.mat');
            end

        elseif w == 2
            if di == 1
                load('mat_files/figureS4B_1.mat');
            elseif di ==2
                load('mat_files/figureS4B_2.mat');
            elseif di == 3
                load('mat_files/figureS4B_3.mat');
            else
                load('mat_files/figureS4B_4.mat');
            end

        else
            if di == 1
                load('mat_files/figureS6F_1.mat');
            elseif di ==2
                load('mat_files/figureS6F_2.mat');
            elseif di == 3
                load('mat_files/figureS6F_3.mat');
            else
                load('mat_files/figureS6F_4.mat');
            end
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
    legend('1/\theta = 0 years','1/\theta = 6 months','1/\theta = 5 years','1/\theta = lifelong','Location','southeast')
    
    if w == 1
        savefig('figure_files/figure2F.fig')
        saveas(gcf,'figure_files/figure2F','epsc')
    elseif w ==2
        savefig('figure_files/figureS4B.fig')
        saveas(gcf,'figure_files/figureS4B','epsc')
    else
        savefig('figure_files/figureS6F.fig')
        saveas(gcf,'figure_files/figureS6F','epsc')
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
