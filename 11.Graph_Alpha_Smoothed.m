%--------------------------------------------------------------------------
% Script to generate Figure A3 (Observed alpha vs. Smoothed alpha) 
% below you can choose between "alpha" and "a" coefficients 
%--------------------------------------------------------------------------

clear; 
clc;
close all

set(0,'defaultTextInterpreter','latex');

% load data

nfVARSpec = '10tc';
nModSpec  = '1';
nMCMCSpec = '1';
modName   = 'SS';  % VAR or SS
sType     = '_factor';    % either empty or _factor

K         = 10;
nKSpec    = ['K',num2str(K),'_'];

% MLE
sName   = ['fVAR',nfVARSpec];
loaddir = [pwd, '/', 'results' ,'/', sName,'/'];
PhatDensCoef = csvread( [loaddir, nKSpec, sName, '_PhatDensCoef', sType, '.csv'], 2, 0); % read table starting two rows below

% Smoothed
sName_s   = ['fVAR',nfVARSpec,'_',modName,nModSpec,'_','MCMC',nMCMCSpec];
loaddir_s = [pwd, '/', 'results' ,'/', sName_s,'/'];
PhatDensCoef_Smoothed = csvread( [loaddir_s, sName_s, '_PhatDensCoef', sType, '_Smoothed.csv'], 1, 0); % read table starting one row below

[T, ~]= size(PhatDensCoef_Smoothed);
period = linspace(1989.25, 1989+0.25*T, T);

figsaveDir = [pwd, '/', 'Figures' ,'/', sName_s,'/'];
[~, ~, ~] = mkdir(figsaveDir);

for kk = 1:K

    figure(kk);clf;
    set(figure(kk),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
    plot(period,PhatDensCoef_Smoothed(:,kk),'Color','r','LineStyle','-','LineWidth',4)
    hold on
    plot(period,PhatDensCoef(:,kk),'Color','b','LineStyle','--','LineWidth',4)
    hold on
    %plot(period, zeros(size(1:T)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
    xlim([period(1) period(end)])
    set(gca,'FontSize',30)
    sNameFile = [sName_s,'_DensCoefs', sType,'_',num2str(kk),'.pdf'];    
    saveas(figure(kk), [figsaveDir sNameFile] );

end
