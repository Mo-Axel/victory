%--------------------------------------------------------------------------
% Script to generate impulse responses
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Housekeeping, Load IRFs of Agg Variables.
%--------------------------------------------------------------------------

clear; 
clc;
close all

set(0,'defaultTextInterpreter','latex');

% set specs 
xmin = 0;
xmax = 3;
xn = 301;
xgrid = linspace(xmin, xmax, xn);

% set specs 
nVARSpec  = '2';
nModSpec  = '1';
nMCMCSpec = '1';
nMod      = 'AltVAR';
n_drawsread = 25; % set how many posterior draws to consider

if nVARSpec == '1'
    nAlt = 'Pctl';
else
    nAlt = 'Gini';
end

sName    = ['AltVAR', nVARSpec, '_MCMC', nMCMCSpec, '_', nAlt];
figDir   = [pwd, '/', 'Figures' ,'/', sName, '/'];
[~, ~, ~] = mkdir(figDir);

% load IRFs

sh_id = 1; % sh_id = 1: TFP shock, sh_id = 2: GDP shock, sh_id = 3: Employment shock

irfDir = [pwd, '/', 'Results' ,'/', sName, '/'];
YY_IRF = csvread( [irfDir, sName, '_IRF_YY_Aggsh',num2str(sh_id),'_pmean.csv'], 1, 0); 
[H, n_all] = size(YY_IRF);
H = H-1;

n_agg = 3;
n_cross = n_all-n_agg;

YY_IRF_uncertainty = zeros(H+1, n_all, n_drawsread);

for pp = 1:n_drawsread
    YY_IRF_uncertainty(:,:,pp)=csvread( [irfDir, sName, '_IRF_YY_AggSh1_', num2str(pp), '.csv'], 1, 0);
end

% first column is TFP growth
TFP_IRF = squeeze(YY_IRF_uncertainty(:,1,:));
TFP_IRF = 100*cumsum(TFP_IRF/400,1);

% second column is GDP growth
GDP_IRF = squeeze(YY_IRF_uncertainty(:,2,:));
GDP_IRF = 100*cumsum(GDP_IRF/400,1);

% third column is employment rate
EMP_IRF = -100*squeeze(YY_IRF_uncertainty(:,3,:));

% last column(s) are summary stats for distribution
Distr_IRF = YY_IRF_uncertainty(:,4:end,:);


%%
%--------------------------------------------------------------------------
% Figure 12: IRFs to a TFP shock
%--------------------------------------------------------------------------

%TFP_IRF_q10 = quantile(TFP_IRF,0.1,2);
%TFP_IRF_q50 = quantile(TFP_IRF,0.5,2);
%TFP_IRF_q90 = quantile(TFP_IRF,0.9,2);

GDP_IRF_q10 = quantile(GDP_IRF,0.1,2);
GDP_IRF_q50 = quantile(GDP_IRF,0.5,2);
GDP_IRF_q90 = quantile(GDP_IRF,0.9,2);

EMP_IRF_q10 = quantile(EMP_IRF,0.1,2);
EMP_IRF_q50 = quantile(EMP_IRF,0.5,2);
EMP_IRF_q90 = quantile(EMP_IRF,0.9,2);

% Figure 9-(i,ii): response of TFP, EMP variables.

hh = 19;

figure(1);clf;
set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(-1:hh,GDP_IRF_q10(1:hh+2,1),'Color','b','LineStyle',':','LineWidth',2)  
hold on
plot(-1:hh,GDP_IRF_q50(1:hh+2,1),'Color','b','LineStyle','-','LineWidth',2) 
hold on
plot(-1:hh,GDP_IRF_q90(1:hh+2,1),'Color','b','LineStyle',':','LineWidth',2) 
set(gca,'FontSize',50)
xlim([-1 hh])

sFigName = [sName, '_IRF_AggSh1_GDP.pdf'];
saveas(figure(1), [figDir sFigName] );
close all

figure(2);clf;
set(figure(2),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(-1:hh,EMP_IRF_q10(1:hh+2,1),'Color','b','LineStyle',':','LineWidth',2)  
hold on
plot(-1:hh,EMP_IRF_q50(1:hh+2,1),'Color','b','LineStyle','-','LineWidth',2) 
hold on
plot(-1:hh,EMP_IRF_q90(1:hh+2,1),'Color','b','LineStyle',':','LineWidth',2) 
set(gca,'FontSize',50)
xlim([-1 hh])

sFigName = [sName, '_IRF_AggSh1_EMP.pdf'];
saveas(figure(2), [figDir '\' sFigName] );
close all


%--------------------------------------------------------------------------
% Figure 12: IRFs to a TFP shock - continued
%--------------------------------------------------------------------------

for ii = 1:n_cross
    
    Pctl_IRF_q10 = quantile(squeeze(Distr_IRF(:,ii,:)), 0.1, 2);
    Pctl_IRF_q50 = quantile(squeeze(Distr_IRF(:,ii,:)), 0.5, 2);
    Pctl_IRF_q90 = quantile(squeeze(Distr_IRF(:,ii,:)), 0.9, 2);
    
    figure(2+ii);clf;
    set(figure(2+ii),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
    plot(-1:hh,Pctl_IRF_q10(1:hh+2,1),'Color','b','LineStyle',':','LineWidth',2)  
    hold on
    plot(-1:hh,Pctl_IRF_q50(1:hh+2,1),'Color','b','LineStyle','-','LineWidth',2) 
    hold on
    plot(-1:hh,Pctl_IRF_q90(1:hh+2,1),'Color','b','LineStyle',':','LineWidth',2)    
    hold on
    plot(-1:hh, Distr_IRF(1,ii,1)*ones(size(-1:hh)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
    set(gca,'FontSize',50)
    xlim([-1 hh])

    if nVARSpec == '1'
     sFigName = [sName, '_IRF_AggSh1_Pctl',num2str(ii),'.pdf'];
    else
     sFigName = [sName, '_IRF_AggSh1_Gini',num2str(ii),'.pdf'];
    end
    
    saveas(figure(2+ii), [figDir '\' sFigName] );
    close all
    
end

