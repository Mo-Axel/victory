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

nfVARSpec = '10tc';
nModSpec  = '1';
nMCMCSpec = '1';
nMod      = 'SS';
n_drawsread = 1000; % set how many posterior draws to consider

sName    = ['fVAR', nfVARSpec, '_', nMod, nModSpec, '_MCMC', nMCMCSpec];
irfDir   = [pwd, '/', 'Results' ,'/', sName, '/'];
figDir   = [pwd, '/', 'Figures' ,'/', sName, '/'];

% load Agg IRFs
YY_IRF = csvread( [irfDir, sName, '_IRF_YY_DistrSh_pmean.csv'], 1, 0); % fVAR IRF
[H, n_all] = size(YY_IRF);
H = H-1;

n_agg = 3;
n_cross = n_all-n_agg;

YY_IRF_uncertainty = zeros(H+1, n_all, n_drawsread);

for pp = 1:n_drawsread
    YY_IRF_uncertainty(:,:,pp)=csvread( [irfDir, sName, '_IRF_YY_DistrSh_', num2str(pp), '.csv'], 1, 0);
end

% first column is TFP growth
TFP_IRF = squeeze(YY_IRF_uncertainty(:,1,:));
TFP_IRF = 100*cumsum(TFP_IRF/400,1);

% second column is GDP growth
GDP_IRF = squeeze(YY_IRF_uncertainty(:,2,:));
GDP_IRF = 100*cumsum(GDP_IRF/400,1);

% third column is employment rate
EMP_IRF = -100*squeeze(YY_IRF_uncertainty(:,3,:));

%%
%--------------------------------------------------------------------------
% Figure 9: IRFs to a TFP shock
%--------------------------------------------------------------------------

TFP_IRF_q10 = quantile(TFP_IRF,0.1,2);
TFP_IRF_q50 = quantile(TFP_IRF,0.5,2);
TFP_IRF_q90 = quantile(TFP_IRF,0.9,2);

GDP_IRF_q10 = quantile(GDP_IRF,0.1,2);
GDP_IRF_q50 = quantile(GDP_IRF,0.5,2);
GDP_IRF_q90 = quantile(GDP_IRF,0.9,2);

EMP_IRF_q10 = quantile(EMP_IRF,0.1,2);
EMP_IRF_q50 = quantile(EMP_IRF,0.5,2);
EMP_IRF_q90 = quantile(EMP_IRF,0.9,2);

% Figure 9-(i,ii): response of TFP, EMP, GDP variables.

hh = 19;

figure(1);clf;
set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(-1:hh,TFP_IRF_q10(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4)  
hold on
plot(-1:hh,TFP_IRF_q50(1:hh+2,1),'Color','b','LineStyle','-','LineWidth',4) 
hold on
plot(-1:hh,TFP_IRF_q90(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4) 
hold on
plot(-1:hh, zeros(hh+2), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
set(gca,'FontSize',50)
xlim([-1 hh])

sFigName = [sName, '_IRF_DistrSh_TFP.pdf'];
saveas(figure(1), [figDir '\' sFigName] );
close all

figure(2);clf;
set(figure(2),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(-1:hh,EMP_IRF_q10(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4)  
hold on
plot(-1:hh,EMP_IRF_q50(1:hh+2,1),'Color','b','LineStyle','-','LineWidth',4) 
hold on
plot(-1:hh,EMP_IRF_q90(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4) 
hold on
plot(-1:hh, zeros(hh+2), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
set(gca,'FontSize',50)
xlim([-1 hh])
ylim([-0.1 0.3])

sFigName = [sName, '_IRF_DistrSh_EMP.pdf'];
saveas(figure(2), [figDir '\' sFigName] );
close all

figure(3);clf;
set(figure(3),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(-1:hh,GDP_IRF_q10(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4)  
hold on
plot(-1:hh,GDP_IRF_q50(1:hh+2,1),'Color','b','LineStyle','-','LineWidth',4) 
hold on
plot(-1:hh,GDP_IRF_q90(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4) 
hold on
plot(-1:hh, zeros(hh+2), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
set(gca,'FontSize',50)
xlim([-1 hh])
ylim([-0.1 0.3])

sFigName = [sName, '_IRF_DistrSh_GDP.pdf'];
saveas(figure(3), [figDir '\' sFigName] );
close all


%%
% Figure 9-(iii-vi): density differential

PhatDens_IRF_uncertainty = zeros(H+1, length(xgrid), n_drawsread);

for pp = 1:n_drawsread
    PhatDens_IRF_uncertainty(:,:,pp)=csvread( [irfDir, sName, '_IRF_PhatDens_DistrSh_', num2str(pp), '.csv'], 1, 0);
end

theta_sinh = 1.0;
ygrid = 1/(2*theta_sinh)*(exp(theta_sinh*xgrid) - exp(-theta_sinh*xgrid));
Jacobian = 1/2*(exp(theta_sinh*xgrid) + exp(-theta_sinh*xgrid));

PhatDens_IRF_ss = PhatDens_IRF_uncertainty(1,:,1)./Jacobian;

PhatDens_IRF_q10 = zeros(H+1,xn);
PhatDens_IRF_q50 = zeros(H+1,xn);
PhatDens_IRF_q90 = zeros(H+1,xn);

for hh=2:(H+1)
    PhatDens_IRF_q10(hh,:)=quantile(squeeze(PhatDens_IRF_uncertainty(hh,:,:)./Jacobian), 0.1, 2); 
    PhatDens_IRF_q50(hh,:)=quantile(squeeze(PhatDens_IRF_uncertainty(hh,:,:)./Jacobian), 0.5, 2); 
    PhatDens_IRF_q90(hh,:)=quantile(squeeze(PhatDens_IRF_uncertainty(hh,:,:)./Jacobian), 0.9, 2); 
end

%PhatDens_IRF = csvread( [irfDir, sName, '_IRF_PhatDens_DistrSh1_pmean.csv'], 1, 0); 

vec_hh = [0; 4];

for ii = 1:length(vec_hh)

hh = vec_hh(ii);

figure(4+ii);clf;
set(figure(4+ii),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(ygrid,PhatDens_IRF_q10(hh+2,:)-PhatDens_IRF_ss,'Color','b','LineStyle','--','LineWidth',4)  
hold on
plot(ygrid,PhatDens_IRF_q50(hh+2,:)-PhatDens_IRF_ss,'Color','b','LineStyle','-','LineWidth',4) 
hold on
plot(ygrid,PhatDens_IRF_q90(hh+2,:)-PhatDens_IRF_ss,'Color','b','LineStyle','--','LineWidth',4) 
hold on
plot(ygrid, zeros(size(ygrid)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)


set(gca,'FontSize',50)
%xlim([ygrid(1) ygrid(end)])
xlim([xmin xmax])

sFigName = [sName, '_IRF_DistrSh_DensDiff_h', num2str(hh), '.pdf'];
saveas(figure(4+ii), [figDir '\' sFigName] );
close all

end


%%
%--------------------------------------------------------------------------
% Figure 10: IRFs to a TFP shock - continued
%--------------------------------------------------------------------------

mean_unrate = 0.051866214912280696; % steady-state value of unemployment rate

cutoff = 1.0;
densvalues_diff = (PhatDens_IRF_ss)./Jacobian;
ygrid_diff      = ygrid(2:xn) - ygrid(1:xn-1);
probmass_diff   = densvalues_diff(2:xn).*ygrid_diff;
baseline = sum(probmass_diff(ygrid(2:xn)<cutoff))+mean_unrate; % add frac zeros too

% Figure 10-(i)
hh = 19;

%BelowCutoff_IRF = csvread( [irfDir, sName, '_IRF_BelowCutoff_DistrSh1_pmean.csv'], 1, 0); 
BelowCutoff_IRF_uncertainty = csvread( [irfDir, sName, '_IRF_BelowCutoff_DistrSh_uncertainty.csv'], 1, 0); 
[~, n_drawsread] = size(BelowCutoff_IRF_uncertainty);

BelowCutoff_IRF_q10 = quantile(BelowCutoff_IRF_uncertainty, 0.1, 2);
BelowCutoff_IRF_q50 = quantile(BelowCutoff_IRF_uncertainty, 0.5, 2);
BelowCutoff_IRF_q90 = quantile(BelowCutoff_IRF_uncertainty, 0.9, 2);

fig_num = 2+length(vec_hh);

figure(fig_num+1);clf;
set(figure(fig_num+1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(-1:hh,baseline+BelowCutoff_IRF_q10(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4)  
hold on
plot(-1:hh,baseline+BelowCutoff_IRF_q50(1:hh+2,1),'Color','b','LineStyle','-','LineWidth',4) 
hold on
plot(-1:hh,baseline+BelowCutoff_IRF_q90(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4) 
hold on
plot(-1:hh, baseline*ones(size(-1:hh)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
set(gca,'FontSize',50)
xlim([-1 hh])

sFigName = [sName, '_IRF_DistrSh_BelowCutoff.pdf'];
saveas(figure(fig_num+1), [figDir '\' sFigName] );
close all


%%
% Figure 10-(ii)
hh = 19;

%Gini_IRF = csvread( [irfDir, sName, '_IRF_Gini_DistrSh1_pmean.csv'], 1, 0); 
Gini_IRF_uncertainty = csvread( [irfDir, sName, '_IRF_Gini_DistrSh_uncertainty.csv'], 1, 0); 
[~, n_drawsread] = size(Gini_IRF_uncertainty);

Gini_IRF_q10 = quantile(Gini_IRF_uncertainty, 0.1, 2);
Gini_IRF_q50 = quantile(Gini_IRF_uncertainty, 0.5, 2);
Gini_IRF_q90 = quantile(Gini_IRF_uncertainty, 0.9, 2);

figure(fig_num+2);clf;
set(figure(fig_num+2),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(-1:hh,Gini_IRF_q10(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4)  
hold on
plot(-1:hh,Gini_IRF_q50(1:hh+2,1),'Color','b','LineStyle','-','LineWidth',4) 
hold on
plot(-1:hh,Gini_IRF_q90(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4)
hold on
plot(-1:hh, Gini_IRF_uncertainty(1,1)*ones(size(-1:hh)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
set(gca,'FontSize',50)
xlim([-1 hh])

sFigName = [sName, '_IRF_DistrSh_Gini.pdf'];
saveas(figure(fig_num+2), [figDir '\' sFigName] );
close all


%%
%--------------------------------------------------------------------------
% Figure 11: IRFs to a TFP shock - continued
%--------------------------------------------------------------------------
hh = 19;

vec_Pctl = [0.1; 0.2; 0.5; 0.8; 0.9];
%Pctl_IRF = csvread( [irfDir, sName, '_IRF_Pctl_DistrSh_pmean.csv'], 1, 0); 
Pctl_IRF_uncertainty = zeros(H+1, length(vec_Pctl), n_drawsread);

for pp = 1:n_drawsread
    Pctl_IRF_uncertainty(:,:,pp)=csvread( [irfDir, sName, '_IRF_Pctl_DistrSh_', num2str(pp), '.csv'], 1, 0);
end


for ii = 1:length(vec_Pctl)
    
    Pctl_IRF_q10 = quantile(squeeze(Pctl_IRF_uncertainty(:,ii,:)), 0.1, 2);
    Pctl_IRF_q50 = quantile(squeeze(Pctl_IRF_uncertainty(:,ii,:)), 0.5, 2);
    Pctl_IRF_q90 = quantile(squeeze(Pctl_IRF_uncertainty(:,ii,:)), 0.9, 2);

    figure(fig_num+2+ii);clf;
    set(figure(fig_num+2+ii),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
    plot(-1:hh,Pctl_IRF_q10(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4)  
    hold on
    plot(-1:hh,Pctl_IRF_q50(1:hh+2,1),'Color','b','LineStyle','-','LineWidth',4) 
    hold on
    plot(-1:hh,Pctl_IRF_q90(1:hh+2,1),'Color','b','LineStyle','--','LineWidth',4)    
    hold on
    plot(-1:hh, Pctl_IRF_uncertainty(1,ii,1)*ones(size(-1:hh)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
    set(gca,'FontSize',50)
    xlim([-1 hh])
    
    if vec_Pctl(ii) == 0.1
        ylim([0.205 0.235])
    elseif vec_Pctl(ii) == 0.9
        ylim([2.4 2.7])
    end        

    sFigName = [sName, '_IRF_DistrSh_Pctl',num2str(ii),'.pdf'];
    saveas(figure(fig_num+2+ii), [figDir '\' sFigName] );
    close all
    
end

