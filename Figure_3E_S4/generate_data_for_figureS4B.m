% Produces mat files that are used in generate_figureS4B.m to generate
% Figure S4B in Chisholm et al., Unravelling the within-host 
% dynamics of Group A Streptococcus from population-level observations of 
% prevalence and strain diversity

clear all
close all

addpath('../Main_code_files')

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool;
end

if ~exist('figure3E_S4_mat_files', 'dir')
    mkdir('figure3E_S4_mat_files')
end

rng(1)

numsims = 16;

R0all = [2 2 2 2];
DIall = [5 7 9 11];

for i = 1:length(R0all)
    R0 = R0all(i);
    DI = DIall(i);
    parfor k = 1: numsims
        [TimeBetweenInfectionSS,~,~,~] = run_model_once_TBI(R0,DI);
        TBI{i,k} = TimeBetweenInfectionSS;
    end
    
    save('figure3E_S4_mat_files/figureS4B.mat','TBI','DIall')
    
end

save('figure3E_S4_mat_files/figureS4B.mat','TBI','DIall')

