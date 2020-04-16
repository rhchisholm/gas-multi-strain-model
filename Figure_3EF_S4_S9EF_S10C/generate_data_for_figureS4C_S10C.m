% Produces mat files that are used in generate_figureS4B.m to generate
% Figure S4C and S10C in Chisholm et al., Unravelling the immune response to Group A 
% Streptococcus infection from population-level observations of 
% prevalence and strain diversity

clear all
close all

addpath('../Main_code_files')

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool;
end

if ~exist('mat_files', 'dir')
    mkdir('mat_files')
end

rng(1)

numsims = 16;

R0all = [2 2 2 2];
DIall = [5 7 9 11];

cross_immunity_effect_on_coinfections = [0 1];

for CIScenario = 1:2
    
    cis = cross_immunity_effect_on_coinfections(CIScenario);
    
    for i = 1:length(R0all)
        
        R0 = R0all(i);
        DI = DIall(i);
        
        parfor k = 1: numsims
            [TimeBetweenInfectionSS,~,~,~] = run_model_once_TBI(R0,DI,cis);
            TBI{i,k} = TimeBetweenInfectionSS;        
        end
        
        if CIScenario == 1
            save('mat_files/figureS4C.mat','TBI','DIall')
        else
            save('mat_files/figureS10C.mat','TBI','DIall')
        end

    end

    if CIScenario == 1
        save('mat_files/figureS4C.mat','TBI','DIall')
    else
        save('mat_files/figureS10C.mat','TBI','DIall')
    end
    
end

