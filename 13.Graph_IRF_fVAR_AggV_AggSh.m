%--------------------------------------------------------------------------
% Script to generate Figure 8 (impulse responses overlay) 
%--------------------------------------------------------------------------

clear; 
clc;
close all

set(0,'defaultTextInterpreter','latex');

% set file/path names
sName1 = 'fVAR10tc_SS1_MCMC1';
sName2 = 'fVAR10tc_SS3_MCMC1';

irfDir1 = [pwd, '/', 'Results' ,'/', sName1, '/'];
irfDir2 = [pwd, '/', 'Results' ,'/', sName2, '/'];
figDir  = [pwd, '/', 'figures' ,'/', 'fVAR10tc_SS1_SS3/'];
[~, ~, ~]  = mkdir(figDir);

% load IRFs
sh_id = 3; % sh_id = 1: TFP shock, sh_id = 2: GDP shock, sh_id = 3: Employment shock
YY_IRF_1 = csvread( [irfDir1, sName1, '_IRF_YY_Aggsh',num2str(sh_id),'_pmean.csv'], 1, 0); 
YY_IRF_2 = csvread( [irfDir2, sName2, '_IRF_YY_Aggsh',num2str(sh_id),'_pmean.csv'], 1, 0); 

[H, n_all1] = size(YY_IRF_1);
[H, n_all2] = size(YY_IRF_2);
H = H-1;
n_drawsread = 1000; % set how many posterior draws to consider

YY_IRF_1_uncertainty = zeros(H+1, n_all1, n_drawsread);
YY_IRF_2_uncertainty = zeros(H+1, n_all2, n_drawsread);

for pp = 1:n_drawsread
    YY_IRF_1_uncertainty(:,:,pp)=csvread( [irfDir1, sName1, '_IRF_YY_AggSh',num2str(sh_id),'_', num2str(pp), '.csv'], 1, 0);
    YY_IRF_2_uncertainty(:,:,pp)=csvread( [irfDir2, sName2, '_IRF_YY_AggSh',num2str(sh_id),'_', num2str(pp), '.csv'], 1, 0);
end


% first column is TFP growth
TFP_IRF_1 = squeeze(YY_IRF_1_uncertainty(:,1,:));
TFP_IRF_1 = 100*cumsum(TFP_IRF_1/400,1);
TFP_IRF_2 = squeeze(YY_IRF_2_uncertainty(:,1,:));
TFP_IRF_2 = 100*cumsum(TFP_IRF_2/400,1);

% second column is GDP growth
GDP_IRF_1 = squeeze(YY_IRF_1_uncertainty(:,2,:));
GDP_IRF_1 = 100*cumsum(GDP_IRF_1/400,1);
GDP_IRF_2 = squeeze(YY_IRF_2_uncertainty(:,2,:));
GDP_IRF_2 = 100*cumsum(GDP_IRF_2/400,1);

% third column is employment rate
EMP_IRF_1 = -100*squeeze(YY_IRF_1_uncertainty(:,3,:));
EMP_IRF_2 = -100*squeeze(YY_IRF_2_uncertainty(:,3,:));

%%
%--------------------------------------------------------------------------
% Figure 9: IRFs to a TFP shock
%--------------------------------------------------------------------------

% first column is TFP growth
TFP_IRF_1_q10 = quantile(TFP_IRF_1,0.1,2);
TFP_IRF_1_q50 = quantile(TFP_IRF_1,0.5,2);
TFP_IRF_1_q90 = quantile(TFP_IRF_1,0.9,2);

TFP_IRF_2_q10 = quantile(TFP_IRF_2,0.1,2);
TFP_IRF_2_q50 = quantile(TFP_IRF_2,0.5,2);
TFP_IRF_2_q90 = quantile(TFP_IRF_2,0.9,2);

% second column is GDP growth
GDP_IRF_1_q10 = quantile(GDP_IRF_1,0.1,2);
GDP_IRF_1_q50 = quantile(GDP_IRF_1,0.5,2);
GDP_IRF_1_q90 = quantile(GDP_IRF_1,0.9,2);

GDP_IRF_2_q10 = quantile(GDP_IRF_2,0.1,2);
GDP_IRF_2_q50 = quantile(GDP_IRF_2,0.5,2);
GDP_IRF_2_q90 = quantile(GDP_IRF_2,0.9,2);

% third column is employment rate
EMP_IRF_1_q10 = quantile(EMP_IRF_1,0.1,2);
EMP_IRF_1_q50 = quantile(EMP_IRF_1,0.5,2);
EMP_IRF_1_q90 = quantile(EMP_IRF_1,0.9,2);

EMP_IRF_2_q10 = quantile(EMP_IRF_2,0.1,2);
EMP_IRF_2_q50 = quantile(EMP_IRF_2,0.5,2);
EMP_IRF_2_q90 = quantile(EMP_IRF_2,0.9,2);


%--------------------------------------------------------------------------
% Plot Responses of AggV to AggSh
%--------------------------------------------------------------------------

figure(1);clf;
set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(1:H+1,TFP_IRF_1_q10,'Color','b','LineStyle','-','LineWidth',4)  
hold on
plot(1:H+1,TFP_IRF_1_q90,'Color','b','LineStyle','-','LineWidth',4) 
hold on
plot(1:H+1,TFP_IRF_2_q10,'Color','r','LineStyle','--','LineWidth',4)  
hold on
plot(1:H+1,TFP_IRF_2_q90,'Color','r','LineStyle','--','LineWidth',4) 
hold on
plot(1:H+1, zeros(size(1:H+1)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
set(gca,'FontSize',50)
xlim([1 (H+1)])

sName = ['Fig_AggSh',num2str(sh_id),'_TFP.pdf'];
saveas(figure(1), [figDir '\' sName] );
close all

figure(2);clf;
set(figure(2),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(1:H+1,GDP_IRF_1_q10,'Color','b','LineStyle','-','LineWidth',4)  
hold on
plot(1:H+1,GDP_IRF_1_q90,'Color','b','LineStyle','-','LineWidth',4) 
hold on
plot(1:H+1,GDP_IRF_2_q10,'Color','r','LineStyle','--','LineWidth',4)  
hold on
plot(1:H+1,GDP_IRF_2_q90,'Color','r','LineStyle','--','LineWidth',4) 
hold on
plot(1:H+1, zeros(size(1:H+1)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
set(gca,'FontSize',50)
xlim([1 (H+1)])

sName = ['Fig_AggSh',num2str(sh_id),'_GDP.pdf'];
saveas(figure(2), [figDir '\' sName] );
close all


figure(3);clf;
set(figure(3),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(1:H+1,EMP_IRF_1_q10,'Color','b','LineStyle','-','LineWidth',4)  
hold on
plot(1:H+1,EMP_IRF_1_q90,'Color','b','LineStyle','-','LineWidth',4) 
hold on
plot(1:H+1,EMP_IRF_2_q10,'Color','r','LineStyle','--','LineWidth',4)  
hold on
plot(1:H+1,EMP_IRF_2_q90,'Color','r','LineStyle','--','LineWidth',4) 
hold on
plot(1:H+1, zeros(size(1:H+1)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
set(gca,'FontSize',50)
xlim([1 (H+1)])

sName = ['Fig_AggSh',num2str(sh_id),'_EMP.pdf'];
saveas(figure(3), [figDir '\' sName] );
close all

