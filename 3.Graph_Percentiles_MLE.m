%--------------------------------------------------------------------------
% Script to overlay percentiles
%--------------------------------------------------------------------------

clear; 
clc;
close all

set(0,'defaultTextInterpreter','latex');

% set specs 
xmin = 0;
xmax = 1;
xn = 101;
xgrid = linspace(xmin, xmax, xn);

nfVARSpec = '10tc';
nKSpec    = 'K4_';
sName = ['fVAR', nfVARSpec];

% load sample percentiles (from the data)
%dataDir = [pwd, '/', 'Data' ,'/'];
dataDir = ['C:/Users/29017/Desktop/Replication/A.FORMAL/CB-fVAR/Data/']
sample_percs = readmatrix([dataDir, 'percentiles_data.csv']); % sample percentiles
sample_percs = sample_percs(2:end,2:end)
%, 2, 1
% load estimated percentiles 
pwd = "C:\Users\29017\Desktop\Replication"
estDir = [pwd, '/', 'Results' ,'/', sName, '/'];
estimated_percs = readmatrix( "C:/Users/29017/Desktop/Replication/results/fVAR10tc/K22_fVAR10tc_PredPctL_MLE.csv")
%estimated_percs = readmatrix( [estDir, nKSpec, sName, '_PredPctl_MLE.csv']); % estimated percentiles
%, 1, 0
%--------------------------------------------------------------------------
% Percentiles Comparison 
%--------------------------------------------------------------------------
T = length(sample_percs(:,1));
figsaveDir = "C:\Users\29017\Desktop\Replication\Figures\fVAR10tc"
%figsaveDir = [pwd, '/', 'Figures' ,'/', sName,'/'];
[~, ~, ~] = mkdir(figsaveDir);

start_period = 1989.25;
period = linspace(start_period,start_period+0.25*(T-1), T);

blue_colors = [0.4 0.4 1; 0 0 0.9; 0 0 0.8; 0 0 0.9; 0.4 0.4 1];
red_colors  = [1 0 0; 0.9 0 0; 0.8 0 0; 0.9 0 0; 1 0 0];

figure(1);clf;
%set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperSize',[11*2 8.5],'PaperPosition',[0.1 0.1 11*2 8.5]);
set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
for ii=1:5
    plot(period, sample_percs(:,ii), 'Color', blue_colors(ii,:), 'LineStyle', '-', 'LineWidth',4)
    hold on
    plot(period, estimated_percs(:,ii), 'Color', red_colors(ii,:), 'LineStyle', '-.', 'LineWidth',4)
    hold on
end
set(gca,'FontSize',30)
xlim([period(1) period(end)])
ylim([0 3])
%axis([1 T 0 3])
sNameFile = [nKSpec, 'Pctl_MLE_All.pdf'];    
saveas(figure(1), [figsaveDir sNameFile] );
%close all


