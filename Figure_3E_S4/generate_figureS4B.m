% Generates Figure S4B in Chisholm et al., Unravelling the within-host 
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

load('figure3E_S4_mat_files/figureS4B.mat');

for j = 1:4
    tbi = [];
    for i = 1:16
        tbi=[tbi;TBI{j,i}];
    end
    tbi(tbi==0)=[];
    alltbi{j} = tbi;
end

figure;
h1=histogram(alltbi{1},'Normalization','pdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[0.07 0.62 1]);
hold on
h2=histogram(alltbi{2},'Normalization','pdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 0.07 0.62]);
hold on
h3=histogram(alltbi{3},'Normalization','pdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[0.62 0.07 1]);
hold on
h4=histogram(alltbi{4},'Normalization','pdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[0.07 0.62 0.62]);
h1.BinWidth = 1;
h2.BinWidth = 1;
h3.BinWidth = 1;
h4.BinWidth = 1;
xlabel('Time until re-infection by the same strain (years)')
ylabel('Frequency')
legend({'1/\theta = 5 years','1/\theta = 7 years','1/\theta = 9 years','1/\theta = 11 years'})
axis([0 70 0 0.12])

savefig('figure3E_S4_figure_files/figureS4B.fig')
saveas(gcf,'figure3E_S4_figure_files/figureS4B','epsc')
