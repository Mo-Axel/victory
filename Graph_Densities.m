%--------------------------------------------------------------------------
% Script to generate Figure 3 (densities) 
%--------------------------------------------------------------------------

clear; 
clc;
close all

set(0,'defaultTextInterpreter','latex');

% load data

dataDir = [pwd, '/', 'Data' ,'/'];
earnings_data = csvread( [dataDir, 'earnings_detrended_inversesign.csv'],1,1);
earnings_t = earnings_data(:,1);
earnings_detrended = earnings_data(:,2);

% set specs and N
% K_vec = [9 10];
K_vec = [10 22];
xmin = 0;
xmax = 3;
xn = 301;
xgrid = linspace(xmin, xmax, xn);

nfVARSpec = '10tc';
sNameDir = ['fVAR',nfVARSpec];
loadDir = [pwd, '/', 'results' ,'/', sNameDir,'/'];
saveDir = [pwd, '/', 'figures' ,'/', sNameDir,'/'];
[~, ~, ~] = mkdir(saveDir);

% given period tt, draw the estimated density

%periods = [110];
periods = [48 88];
%periods = 1:115;


for pp = 1:length(periods)
    
    tt = periods(pp);
    timeidx = 1989+0.25*(tt-1);

    selecteddraws_t = earnings_detrended(earnings_t==timeidx);

    color_choice = ['r','b','k','m','g'];
    linestyle_choice = ['-',':','-','--','-'];

    figure(1);clf;
    set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
    histogram(selecteddraws_t, 40, 'Normalization', 'pdf','FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
    % we previously used 20 bins, 10 bins gives more smoothing, and 40 bins
    % highlights the peaks that the log spline estimator is picking up
    hold on

    for ii = 1:length(K_vec)

        K =  K_vec(ii);
        sNameFile = ['K',num2str(K_vec(ii)),'_fVAR',nfVARSpec];    
        PhatDensValue = csvread( [loadDir, sNameFile, '_PhatDensValue.csv'], 1, 0); % read table starting one row below

        plot(xgrid,PhatDensValue(tt,:),'Color',color_choice(ii),'LineStyle',linestyle_choice(ii),'LineWidth',4)
        hold on

    end
       
    legend('histogram','K=10','K=22', 'Location','northeast')    
    set(gca,'FontSize',30)
    
    sName = ['fVAR',nfVARSpec,'_Density_Period',num2str(periods(pp)),'.pdf'];
    saveas(figure(1), [saveDir sName] );
    % close all


end




