% Generates heat maps in Figure S1BFJNRV 

close all
clear all

set(0,'DefaultTextFontName','Arial')
set(0,'DefaultTextFontSize',33)
set(0,'DefaultAxesFontSize',33)
set(0,'DefaultAxesFontName','Arial')

if ~exist('figure_files', 'dir')
    mkdir('figure_files')
end

l = ['ABCD';'EFGH';'IJKL';
    'MNOP';'QRST';'UVWX'];

set(0, 'DefaultFigureRenderer', 'painters');

for w = 1 : 6
    
    load(sprintf('mat_files/figureS1_%s.mat',l(w,:)))

    
    [xm,ym] = size(SSPrev);
    x = (1:xm)-1;
    y = (1:ym)/365;
    SSPrev(SSPrev == 0) = NaN;

    figure
    plotheatmap(y,x,SSPrev')
    axis([8 10 0.4 42.6])
    
    savefig(sprintf('figure_files/figureS1_%s.fig',l(w,2)))
    saveas(gcf,sprintf('figure_files/figureS1_%s',l(w,2)),'epsc')
    
end

%savefig('figure_files/figure1B.fig')
%saveas(gcf,'figure_files/figure1B','epsc')

function D = div(SSP)

    SSP1 = SSP - 1;
    N = sum(SSP);
    D = N .* (N - 1);
    sumSSP = sum(SSP .* SSP1);
    D = D ./ sumSSP;
    D(D == Inf) = N(D == Inf); 
    D(isnan(D)) = 0;

end

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

    pcolor(x,y,zz');
    shading flat
    h = colorbar;
    %ylabel(h, {'Number of observed infections'})
    caxis([0 500])
    xlabel('Time (years)')
    ylabel('Strain number')
    % create custom color map
    map = linspace(0,500,100);
    map = fliplr(map)';
    map = [map map map]/500;
    colormap(map)

end