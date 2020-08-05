% Generates Figure 2A-D, S2, S3, S6A-D and S7A-BB

clear all
close all

set(0,'DefaultTextFontName','Arial')
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultAxesFontName','Arial')

if ~exist('figure_files', 'dir')
    mkdir('figure_files')
end

% Scenario w = 1:2 : figure2A-D
% Scenario w = 3:4 : figure S2
% Scenario w = 5:6 : figure S6A-D


for w = 1 : 6
    
    if w == 1
        load('mat_files/figure2_AB.mat');
    elseif w == 2
        load('mat_files/figure2_CD.mat');
    elseif w == 3       
        load('mat_files/figureS2_AB.mat');
    elseif w == 4
        load('mat_files/figureS2_CD.mat');
    elseif w == 5
        load('mat_files/figureS6_AB.mat');
    else
        load('mat_files/figureS6_CD.mat');
    end
    
    for i = 1 : length(StrainSpecific)

        for j = 1 : length(CrossStrain)

            SSPrevTimePoint(:,:,i,j) = TimePrev(j,i,:,:);
            NumAgentsInfected(:,:,i,j) = TimeAgentsInfectedByKStrains(j,i,:,:);

            for k = 1 : NumberSimulations

                TotalPrev(i,j,k) = sum(SSPrevTimePoint(:,k,i,j),1);
                TotalAgentsInfected(i,j,k) = sum(NumAgentsInfected(:,k,i,j),1);

            end

        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE DIVERSITY, PREVALENC, COINFECTION %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    % Calculate final diversity of each simulation
    Diversity = squeeze(div(SSPrevTimePoint));
    
    % Calculate mean final diversity
    z = mean(Diversity,1);
    z = squeeze(z);

    % Calculate mean final prevalence of each simulation
    AvgTotalAgents = mean(TotalAgentsInfected,3);
    Prevalence = AvgTotalAgents / Nagents;
     
    % Grid of strengths of strain-specific and cross-strain immunity
    x = StrainSpecific;
    dx = (x(2) - x(1))/2;
    y = CrossStrain;
    dy = (y(2) - y(1))/2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MEAN FINAL DIVERSITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Color is white where diversity does not exist (or where parameter 
    % space is not sampled)   
    z(z==0) = NaN; 
    
    % Otherwise color reflects diversity.  Need to adjust boundary
    zz = zeros(size(z,1)+1,size(z,2)+1);
    zz(1:size(z,1),1:size(z,2))=z;
    zz(1:size(z,1),end)=z(:,end);
    zz(end,1:size(z,2))=z(end,:);
    
    % Reposition data points so that the middle of a pixel aligns value of 
    % s-s or c-s immunity strength
    figure
    pcolor([y-dy y(end)+dy],[x-dx x(end)+dx],zz);
    shading flat
    colorbar
    xlabel('Strength of strain-specific immunity')
    ylabel('Strength of cross-strain immunity')
    if w <= 4
        caxis([0 32.6153])
    else
        caxis([0 32.20525])
    end
    
    if w == 1
        savefig('figure_files/figure2_A.fig')
        saveas(gcf,'figure_files/figure2_A','epsc')
    elseif w == 2
        savefig('figure_files/figure2_C.fig')
        saveas(gcf,'figure_files/figure2_C','epsc')
    elseif w == 3
        savefig('figure_files/figureS2_A.fig')
        saveas(gcf,'figure_files/figureS2_A','epsc')
    elseif w == 4
        savefig('figure_files/figureS2_C.fig')
        saveas(gcf,'figure_files/figureS2_C','epsc')
    elseif w == 5
        savefig('figure_files/figureS6_A.fig')
        saveas(gcf,'figure_files/figureS6_A','epsc')
    else
        savefig('figure_files/figureS6_C.fig')
        saveas(gcf,'figure_files/figureS6_C','epsc')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MEAN FINAL PREVALENCE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Color is white where parameter space is not sampled
    for i = 1 : length(StrainSpecific)
        for j = 1 : i - 1      
               Prevalence(i,j) = NaN;  
        end
    end
    
    % Color reflects prevalence.  Need to adjust boundary
    zz = zeros(size(z,1)+1,size(z,2)+1);
    zz(1:size(z,1),1:size(z,2))=Prevalence;
    zz(1:size(z,1),end)=Prevalence(:,end);
    zz(end,1:size(z,2))=Prevalence(end,:);
    
    % Reposition data points so that the middle of a pixel aligns value of 
    % s-s or c-s immunity strength
    figure
    pcolor([y-dy y(end)+dy],[x-dx x(end)+dx],zz*100);
    shading flat
    colorbar
    if w <= 4
        caxis([0 71.492])
    else 
        caxis([0 71.664])
    end
    xlabel('Strength of strain-specific immunity')
    ylabel('Strength of cross-strain immunity')
    
    if w == 1
        savefig('figure_files/figure2_B.fig')
        saveas(gcf,'figure_files/figure2_B','epsc')
    elseif w == 2
        savefig('figure_files/figure2_D.fig')
        saveas(gcf,'figure_files/figure2_D','epsc')
    elseif w == 3
        savefig('figure_files/figureS2_B.fig')
        saveas(gcf,'figure_files/figureS2_B','epsc')
    elseif w == 4
        savefig('figure_files/figureS2_D.fig')
        saveas(gcf,'figure_files/figureS2_D','epsc')
    elseif w == 5
        savefig('figure_files/figureS6_B.fig')
        saveas(gcf,'figure_files/figureS6_B','epsc')
    else
        savefig('figure_files/figureS6_D.fig')
        saveas(gcf,'figure_files/figureS6_D','epsc')
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUPPLEMENTARY FIGURE 3: CALCULATE DIFFERENCE IN DIVERSITY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenario w = 1 : figure S3A
% Scenario w = 2 : figure S3B
% Scenario w = 3 : figure S10A
% Scenario w = 4 : figure S10B

for w = 1 : 4
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE DIVERSITY SHORT DURATION IMMUNITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if w == 1
        load('mat_files/figure2_AB.mat');
    elseif w == 2
        load('mat_files/figureS2_AB.mat');
    elseif w == 3
        load('mat_files/figureS6_AB.mat');
    else
        load('mat_files/figureS7_AB.mat');
    end

    for i = 1 : length(StrainSpecific)
        for j = 1 : length(CrossStrain)
            SSPrevTimePoint(:,:,i,j) = TimePrev(j,i,:,:);
        end
    end
     
    Diversity = squeeze(div(SSPrevTimePoint));
    
    % Calculate mean final diversity
    z = mean(Diversity,1);
    z1 = squeeze(z);
    
    % Color is white where diversity does not exist (or where parameter 
    % space is not sampled)   
    z1(z1==0) = NaN; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE DIVERSITY LONG DURATION IMMUNITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if w == 1
        load('mat_files/figure2_CD.mat');
    elseif w == 2
        load('mat_files/figureS2_CD.mat');
    elseif w == 3
        load('mat_files/figureS6_CD.mat');
    else
        load('mat_files/figureS7_CD.mat');
    end
    
    for i = 1 : length(StrainSpecific)
        for j = 1 : length(CrossStrain)
            SSPrevTimePoint(:,:,i,j) = TimePrev(j,i,:,:);
        end
    end

    Diversity = squeeze(div(SSPrevTimePoint));
    
    % Calculate mean final diversity
    z = mean(Diversity,1);
    z2 = squeeze(z);
    
    % Color is white where diversity does not exist (or where parameter 
    % space is not sampled)   
    z2(z2==0) = NaN; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE DIFFERENCE IN DIVERSITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z = z2 - z1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MEAN FINAL DIVERSITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Grid of strengths of strain-specific and cross-strain immunity
    x = StrainSpecific;
    dx = (x(2) - x(1))/2;
    y = CrossStrain;
    dy = (y(2) - y(1))/2;
    
    % Need to adjust boundary
    zz = zeros(size(z,1)+1,size(z,2)+1);
    zz(1:size(z,1),1:size(z,2))=z;
    zz(1:size(z,1),end)=z(:,end);
    zz(end,1:size(z,2))=z(end,:);
    
    % Reposition data points so that the middle of a pixel aligns value of 
    % s-s or c-s immunity strength
    figure
    mymap1 = [ones(20,1) linspace(0,1,20)' linspace(0,1,20)'];
    mymap2 = [flip(linspace(0,1,75)') flip(linspace(0,1,75)') ones(75,1)];
    mymap = [mymap1;mymap2];
    
    s=pcolor([y-dy y(end)+dy],[x-dx x(end)+dx],zz);
    caxis([-2 7.5])
    colormap(mymap);
    shading flat
    colorbar
    xlabel('Strength of strain-specific immunity')
    ylabel('Strength of cross-strain immunity')
  
    if w == 1
        savefig('figure_files/figureS3_A.fig')
        saveas(gcf,'figure_files/figureS3_A','epsc')
    elseif w == 2
        savefig('figure_files/figureS3_B.fig')
        saveas(gcf,'figure_files/figureS3_B','epsc')
    elseif w == 3
        savefig('figure_files/figureS7_A.fig')
        saveas(gcf,'figure_files/figureS7_A','epsc')
    else
        savefig('figure_files/figureS7_B.fig')
        saveas(gcf,'figure_files/figureS7_B','epsc')
    end

end

% Function that calucations period diversity from matrix of number of
% infections of each strain at each time step (SSP)
function D = div(SSP)

    SSP1 = SSP - 1;
    N = sum(SSP);
    D = N .* (N - 1);
    sumSSP = sum(SSP .* SSP1);
    D = D ./ sumSSP;
    D(D == Inf) = N(D == Inf);
    D(isnan(D)) = 0;

end