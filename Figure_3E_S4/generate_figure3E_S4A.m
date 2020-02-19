% Generates Figure 3E and S4A in Chisholm et al., Unravelling the within-host 
% dynamics of Group A Streptococcus from population-level observations of 
% prevalence and strain diversity

clear all
close all

set(0,'DefaultTextFontName','Arial')
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultAxesFontName','Arial')

if ~exist('figure3E_S4_figure_files', 'dir')
    mkdir('figure3E_S4_figure_files')
end

for w = 1 : 2
    
    %%%%%%%%%%%%%
    % LOAD DATA %
    %%%%%%%%%%%%%
    if w == 1
        load('figure3E_S4_mat_files/figure3E.mat');
    else
        load('figure3E_S4_mat_files/figureS4A.mat');
    end
    
    R0 = 0.0301 / (1/71/52.14 + 1/2 + params(9) / params(8));
    
    TotalAgentsInfected = zeros(1,length(Contacts),NumberSimulations);
    TotalPrev = zeros(1,length(Contacts),NumberSimulations);
    Diversity = zeros(1,length(Contacts),NumberSimulations);
    
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
    
    if w == 1
        %title('Low resistance to co-infection')
        legend('1/\theta = 0 years','1/\theta = 6 months','1/\theta = 5 years','1/\theta = lifelong','Location','southeast')
        savefig('figure3E_S4_figure_files/figure3E.fig')
        saveas(gcf,'figure3E_S4_figure_files/figure3E','epsc')
    else 
        %title('High resistance to co-infection')
        legend('1/\theta = 0 years','1/\theta = 6 months','1/\theta = 5 years','1/\theta = lifelong','Location','southeast')
        savefig('figure3E_S4_figure_files/figureS4A.fig')
        saveas(gcf,'figure3E_S4_figure_files/figureS4A','epsc')
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
