% Generates Figure S1 in Chisholm et al., Unravelling the immune response 
% to Group A Streptococcus infection from population-level observations of 
% prevalence and strain diversity

close all
clear all

set(0,'DefaultTextFontName','Arial')
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultAxesFontName','Arial')

if ~exist('figure_files', 'dir')
    mkdir('figure_files')
end

l = ['ABC';'DEF';'GHI'];

set(0, 'DefaultFigureRenderer', 'painters');

for w = 1 : 3
    
    load(sprintf('mat_files/figureS1_%s.mat',l(w,:)))
    
    time = linspace(0,DurationSimulation,ceil(DurationSimulation*52.14*7));

   % Plot number of infections of each strain as function of time
    figure
    xlabel('Time (years)')
    line(time,SSPrev)
    axis([0 time(end) 0 500])
    ylabel('Number of infections')
    
    savefig(sprintf('figure_files/figureS1_%s.fig',l(w,1)))
    saveas(gcf,sprintf('figure_files/figureS1_%s',l(w,1)),'epsc')
    
    % Plot number of infections per host at end of simulation
    figure
    CoinfDistribution = AgentsInfectedByKStrains(1:6,end);
    total5 = sum(CoinfDistribution);
    bar(1:6,CoinfDistribution/total5)
    xlabel('Number of infections per host')
    ylabel('Final host proportion')
    axis([0 7 0 1])
    
    savefig(sprintf('figure_files/figureS1_%s.fig',l(w,2)))
    saveas(gcf,sprintf('figure_files/figureS1_%s',l(w,2)),'epsc')
 
    % Plot diversity through time
    Diversity = div(SSPrev);
    %Diversity rolling average over 10 years from prevalence data
    for i = 1:length(time)
        timesteps = 1:1:i;
        PrevD(:,i)=sum(SSPrev(:,timesteps),2);
    end

    Diversity2 = div(PrevD);

    Prevalence = sum(AgentsInfectedByKStrains,1)/Nagents;
        
    figure
    line(time,Prevalence*100,'Color','k')
    ylabel('Prevalence')
    axis([0 time(end) 0 100])
    xlabel('Time (years)')
    yyaxis right
    line(time,Diversity2,'LineWidth',2,'Color','b')
    axis([0 time(end) 0 45])
    ylabel('Diversity')
    
    savefig(sprintf('figure_files/figureS1_%s.fig',l(w,3)))
    saveas(gcf,sprintf('figure_files/figureS1_%s',l(w,3)),'epsc')
    

end

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