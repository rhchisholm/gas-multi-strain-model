% Run model once

clear all
close all

set(0,'DefaultTextFontName','Arial')
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultAxesFontName','Arial')

rng(3) % seed for the random numbers, needed to reproduce a simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Duration of simulation (years)
DurationSimulation = 20;
% Number of strains
Nstrains = 42;
% Total number of hosts
Nagents = 2500; 
% Age that all hosts die (years)
AgeDeath = 71; 
% Basic Reproduction Number 
BasicReproductionNumber = 2.07;  
% Duration of immunity (weeks)
Dimmunity = 0.5 * 52.14; 
% Resistance to co-infection
x = 10;
% Migration rate per week per population
alpha = 3;  % daily per capita rate is alpha / 7 / Nagents
% Number of contacts per week
Cperweek = 34.53;
% Strength of strain-specific immunity
sigma = 1;
% Strength of cross-strain immunity
omega = 0.1;

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

[SSPrev,AgentsInfectedByKStrains] = ...
                simulator(AgentCharacteristics, ImmuneStatus, params, 0, 1);
           
time = linspace(0,DurationSimulation,size(SSPrev,2));

% Plot number of infections of each strain as function of time
figure
xlabel('Time (years)')
line(time,SSPrev,'LineWidth',2)
%axis([DurationSimulation-2 DurationSimulation 0 500])
ylabel('Number of infections')


% Plot number of infections per host at end of simulation
figure
CoinfDistribution = AgentsInfectedByKStrains(1:6,end);
total5 = sum(CoinfDistribution);
bar(1:6,CoinfDistribution/total5)
xlabel('Number of infections per host')
ylabel('Final host proportion')
axis([0 7 0 1])


% Plot diversity through time
Diversity = div(SSPrev);
Prevalence = sum(AgentsInfectedByKStrains,1)/Nagents;

figure
line(time,Prevalence*100,'Color','k','LineWidth',2)
ylabel('Prevalence')
axis([DurationSimulation-2 DurationSimulation 0 20])
xlabel('Time (years)')
yyaxis right
line(time,Diversity,'LineWidth',2)
axis([DurationSimulation-2 DurationSimulation 0 10])
ylabel('Diversity')


Outbreak = SSPrev;       
[xm,ym] = size(Outbreak);
xx = (1:xm);
y = (1:ym)/365;
figure
plotheatmap(xx,y,Outbreak)
%axis([DurationSimulation-2 DurationSimulation 0.4 42.6])
xlabel('Time (years)')
ylabel('Strain number')

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
figure
plot(time,DivAvg2,'LineWidth',2)
ylabel('10-year rolling avg diversity (prev)')


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

% Plots heat map showing extant strains
function plotheatmap(x,y,z)

    z(z == 0) = NaN;

    % Grid of strain numbers and months
    dx = (x(2) - x(1))/2;
    dy = (y(2) - y(1))/2;

    % Reposition data points so that the middle of a pixel aligns value  
    % of x and y
    x = [x-dx x(end)+dx];
    y = [y-dy y(end)+dy];

    % Adjust boundary
    zz = zeros(size(z,1)+1,size(z,2)+1);
    zz(1:size(z,1),1:size(z,2))=z;
    zz(1:size(z,1),end)=z(:,end);
    zz(end,1:size(z,2))=z(end,:);

    pcolor(y,x,zz);
    shading flat
    colorbar
    caxis([0 90])
    % create custom color map
    map = linspace(0,1,100);
    map = fliplr(map)';
    map = [map map map];
    colormap(map)
    axis([-0.03 10.005 0.4 42.6])

end






