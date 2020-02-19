% Generates Figure 3A-D, S2 and S3 in Chisholm et al., Unravelling the 
% within-host dynamics of Group A Streptococcus from population-level 
% observations of prevalence and strain diversity

clear all
close all

set(0,'DefaultTextFontName','Arial')
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultAxesFontName','Arial')

if ~exist('figure3ABCD_S2_S3_figure_files', 'dir')
    mkdir('figure3ABCD_S2_S3_figure_files')
end

for w = 1 : 4
    
    if w == 1
        load('figure3ABCD_S2_S3_mat_files/figure3_AB.mat');
    elseif w == 2
        load('figure3ABCD_S2_S3_mat_files/figure3_CD.mat');
    elseif w == 3       
        load('figure3ABCD_S2_S3_mat_files/figureS2_AB.mat');
    else
        load('figure3ABCD_S2_S3_mat_files/figureS2_CD.mat');
    end
    
    for i = 1 : length(StrainSpecific)

        for j = 1 : length(CrossStrain)

            SSPrevTimePoint(:,:,i,j) = TimePrev(j,i,:,:);
            NumAgentsInfected(:,:,i,j) = TimeAgentsInfectedByKStrains(j,i,:,:);

            for k = 1 : NumberSimulations

                %TotalPrev(i,j,k) = sum(SSPrevTimePoint(:,k,i,j),1);
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
    caxis([0 32.6153])
    
    if w == 1
        savefig('figure3ABCD_S2_S3_figure_files/figure3_A.fig')
        saveas(gcf,'figure3ABCD_S2_S3_figure_files/figure3_A','epsc')
    elseif w == 2
        savefig('figure3ABCD_S2_S3_figure_files/figure3_C.fig')
        saveas(gcf,'figure3ABCD_S2_S3_figure_files/figure3_C','epsc')
    elseif w == 3
        savefig('figure3ABCD_S2_S3_figure_files/figureS2_A.fig')
        saveas(gcf,'figure3ABCD_S2_S3_figure_files/figureS2_A','epsc')
    else
        savefig('figure3ABCD_S2_S3_figure_files/figureS2_C.fig')
        saveas(gcf,'figure3ABCD_S2_S3_figure_files/figureS2_C','epsc')
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
    caxis([0 71.492])
    xlabel('Strength of strain-specific immunity')
    ylabel('Strength of cross-strain immunity')
    
    if w == 1
        savefig('figure3ABCD_S2_S3_figure_files/figure3_B.fig')
        saveas(gcf,'figure3ABCD_S2_S3_figure_files/figure3_B','epsc')
    elseif w == 2
        savefig('figure3ABCD_S2_S3_figure_files/figure3_D.fig')
        saveas(gcf,'figure3ABCD_S2_S3_figure_files/figure3_D','epsc')
    elseif w == 3
        savefig('figure3ABCD_S2_S3_figure_files/figureS2_B.fig')
        saveas(gcf,'figure3ABCD_S2_S3_figure_files/figureS2_B','epsc')
    else
        savefig('figure3ABCD_S2_S3_figure_files/figureS2_D.fig')
        saveas(gcf,'figure3ABCD_S2_S3_figure_files/figureS2_D','epsc')
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUPPLEMENTARY FIGURE 2: CALCULATE DIFFERENCE IN DIVERSITY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for w = 1 : 2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE DIVERSITY SHORT DURATION IMMUNITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if w == 1
        load('figure3ABCD_S2_S3_mat_files/figure3_AB.mat');
    elseif w == 2
        load('figure3ABCD_S2_S3_mat_files/figureS2_AB.mat');
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
        load('figure3ABCD_S2_S3_mat_files/figure3_CD.mat');
    elseif w == 2
        load('figure3ABCD_S2_S3_mat_files/figureS2_CD.mat');
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
        savefig('figure3ABCD_S2_S3_figure_files/figureS3_A.fig')
        saveas(gcf,'figure3ABCD_S2_S3_figure_files/figureS3_A','epsc')
    else
        savefig('figure3ABCD_S2_S3_figure_files/figureS3_B.fig')
        saveas(gcf,'figure3ABCD_S2_S3_figure_files/figureS3_B','epsc')
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