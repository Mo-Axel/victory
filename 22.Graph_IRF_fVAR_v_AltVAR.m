%--------------------------------------------------------------------------
% Script to generate fVAR vs. AltVAR (impulse responses overlay) 
%--------------------------------------------------------------------------

clear; 
clc;
close all

set(0,'defaultTextInterpreter','latex');

selectGini = 0;

% set file/path names
sName1 = 'fVAR10tc_SS1_MCMC1';

if selectGini == 1
    sName2 = 'AltVAR2_MCMC1_Gini';
    figDir  = [pwd, '/', 'figures' ,'/', 'fVAR10tc_AltVAR2_Gini/'];    
else
    sName2 = 'AltVAR1_MCMC1_Pctl';
    figDir  = [pwd, '/', 'figures' ,'/', 'fVAR10tc_AltVAR1_Pctl/'];
end
[~, ~, ~]  = mkdir(figDir);


irfDir1 = [pwd, '/', 'Results' ,'/', sName1, '/'];
irfDir2 = [pwd, '/', 'Results' ,'/', sName2, '/'];


% load IRFs
sh_id = 1; 
YY_IRF_1 = csvread( [irfDir1, sName1, '_IRF_YY_Aggsh',num2str(sh_id),'_pmean.csv'], 1, 0); 
YY_IRF_2 = csvread( [irfDir2, sName2, '_IRF_YY_Aggsh',num2str(sh_id),'_pmean.csv'], 1, 0); 

[H_1, n_all_1] = size(YY_IRF_1); %% NOTE: H differs from H_2
[H_2, n_all_2] = size(YY_IRF_2);
H_1 = H_1-1;
H_2 = H_2-1;

n_agg = 3;

n_drawsread = 1000; % set how many posterior draws to consider

YY_IRF_1_uncertainty = zeros(H_1+1, n_all_1, n_drawsread);
YY_IRF_2_uncertainty = zeros(H_2+1, n_all_2, n_drawsread);

for pp = 1:n_drawsread
    YY_IRF_1_uncertainty(:,:,pp)=csvread( [irfDir1, sName1, '_IRF_YY_AggSh',num2str(sh_id),'_', num2str(pp), '.csv'], 1, 0);
    YY_IRF_2_uncertainty(:,:,pp)=csvread( [irfDir2, sName2, '_IRF_YY_AggSh',num2str(sh_id),'_', num2str(pp), '.csv'], 1, 0);
end

H = min(H_1,H_2);

% first column is TFP growth
TFP_IRF_1 = squeeze(YY_IRF_1_uncertainty(1:H+1,1,:));
TFP_IRF_1 = 100*cumsum(TFP_IRF_1/400,1);
TFP_IRF_2 = squeeze(YY_IRF_2_uncertainty(1:H+1,1,:));
TFP_IRF_2 = 100*cumsum(TFP_IRF_2/400,1);

% second column is GDP growth
GDP_IRF_1 = squeeze(YY_IRF_1_uncertainty(1:H+1,2,:));
GDP_IRF_1 = 100*cumsum(GDP_IRF_1/400,1);
GDP_IRF_2 = squeeze(YY_IRF_2_uncertainty(1:H+1,2,:));
GDP_IRF_2 = 100*cumsum(GDP_IRF_2/400,1);

% third column is employment rate
EMP_IRF_1 = -100*squeeze(YY_IRF_1_uncertainty(1:H+1,3,:));
EMP_IRF_2 = -100*squeeze(YY_IRF_2_uncertainty(1:H+1,3,:));

%%
%--------------------------------------------------------------------------
% Figure 13-14: IRFs to a TFP shock
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

%%
%--------------------------------------------------------------------------
% Plot Responses of Distr to AggSh
%--------------------------------------------------------------------------

H = min(H_1,H_2);

if selectGini == 1

    % FVAR
    Gini_IRF_1_uncertainty = csvread( [irfDir1, sName1, '_IRF_Gini_AggSh1_uncertainty.csv'], 1, 0);
    BelowCutoff_IRF_1_uncertainty = csvread( [irfDir1, sName1, '_IRF_BelowCutoff_AggSh1_uncertainty.csv'], 1, 0);
     
    Gini_IRF_1 = Gini_IRF_1_uncertainty(1:H+1,:)
    BelowCutoff_IRF_1 = BelowCutoff_IRF_1_uncertainty(1:H+1,:)

    Gini_IRF_1_q10 = quantile(Gini_IRF_1,0.1,2);
    Gini_IRF_1_q50 = quantile(Gini_IRF_1,0.5,2);
    Gini_IRF_1_q90 = quantile(Gini_IRF_1,0.9,2);

    BelowCutoff_IRF_1_q10 = quantile(BelowCutoff_IRF_1,0.1,2);
    BelowCutoff_IRF_1_q50 = quantile(BelowCutoff_IRF_1,0.5,2);
    BelowCutoff_IRF_1_q90 = quantile(BelowCutoff_IRF_1,0.9,2);
    
    % AltVAR: column 4 is Gini, column 5 is BelowCutoff    
    Gini_IRF_2        = Gini_IRF_1(1,1) + squeeze(YY_IRF_2_uncertainty(1:H+1,4,:));
    BelowCutoff_IRF_2 = squeeze(YY_IRF_2_uncertainty(1:H+1,5,:));

    Gini_IRF_2_q10 = quantile(Gini_IRF_2,0.1,2);
    Gini_IRF_2_q50 = quantile(Gini_IRF_2,0.5,2);
    Gini_IRF_2_q90 = quantile(Gini_IRF_2,0.9,2);

    BelowCutoff_IRF_2_q10 = quantile(BelowCutoff_IRF_2,0.1,2);
    BelowCutoff_IRF_2_q50 = quantile(BelowCutoff_IRF_2,0.5,2);
    BelowCutoff_IRF_2_q90 = quantile(BelowCutoff_IRF_2,0.9,2);
    
    figure(4);clf;
    set(figure(4),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
    plot(1:H+1, Gini_IRF_1_q10, 'Color', 'b', 'LineStyle', '-', 'LineWidth',4)
    hold on
    plot(1:H+1, Gini_IRF_1_q90, 'Color', 'b', 'LineStyle', '-', 'LineWidth',4)
    hold on
    plot(1:H+1, Gini_IRF_2_q10, 'Color', 'r', 'LineStyle', '--', 'LineWidth',4)
    hold on
    plot(1:H+1, Gini_IRF_2_q90, 'Color', 'r', 'LineStyle', '--', 'LineWidth',4)
    hold on    
    plot(1:H+1, Gini_IRF_1(1,1)*ones(size(1:H+1)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
    set(gca,'FontSize',50)
    xlim([1 (H+1)])

    sName = ['Fig_AggSh',num2str(sh_id),'_Gini.pdf'];
    saveas(figure(4), [figDir '\' sName] );
    close all


    figure(5);clf;
    set(figure(5),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
    plot(1:H+1, BelowCutoff_IRF_1_q10, 'Color', 'b', 'LineStyle', '-', 'LineWidth',4)
    hold on
    plot(1:H+1, BelowCutoff_IRF_1_q90, 'Color', 'b', 'LineStyle', '-', 'LineWidth',4)
    hold on
    plot(1:H+1, BelowCutoff_IRF_2_q10, 'Color', 'r', 'LineStyle', '--', 'LineWidth',4)
    hold on
    plot(1:H+1, BelowCutoff_IRF_2_q90, 'Color', 'r', 'LineStyle', '--', 'LineWidth',4)
    hold on
    plot(1:H+1, zeros(size(1:H+1)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
    set(gca,'FontSize',50)
    xlim([1 (H+1)])

    sName = ['Fig_AggSh',num2str(sh_id),'_BelowCutoff.pdf'];
    saveas(figure(5), [figDir '\' sName] );
    close all

else
    
    % fVAR
    vec_Pctl = [0.1; 0.2; 0.5; 0.8; 0.9];
    n_percs  = length(vec_Pctl);
    
    Pctl_IRF1_uncertainty = zeros(H+1, n_percs, n_drawsread);

    for pp = 1:n_drawsread
        Pctl_IRF1_uncertainty(:,:,pp)=csvread( [irfDir1, sName1, '_IRF_Pctl_AggSh',num2str(sh_id),'_', num2str(pp), '.csv'], 1, 0);
    end
    
    for ii = 1:n_percs

        set(figure(n_agg+ii),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);

        Pctl_IRF1_q10 = quantile(squeeze(Pctl_IRF1_uncertainty(:,ii,:)), 0.1, 2);
        Pctl_IRF1_q90 = quantile(squeeze(Pctl_IRF1_uncertainty(:,ii,:)), 0.9, 2);
 
        plot(1:H+1, Pctl_IRF1_q10, 'Color', 'b', 'LineStyle', '-', 'LineWidth',4)
        hold on
        plot(1:H+1, Pctl_IRF1_q90, 'Color', 'b', 'LineStyle', '-', 'LineWidth',4)
        hold on

        Pctl_IRF2     = squeeze(YY_IRF_2_uncertainty(1:H+1,n_agg+ii,:));
        Pctl_IRF2_q10 = quantile(Pctl_IRF2, 0.1, 2);
        Pctl_IRF2_q90 = quantile(Pctl_IRF2, 0.9, 2);

        plot(1:H+1, (Pctl_IRF1_q10(1)+Pctl_IRF2_q10), 'Color', 'r', 'LineStyle', '--', 'LineWidth',4)
        hold on
        plot(1:H+1, (Pctl_IRF1_q10(1)+Pctl_IRF2_q90), 'Color', 'r', 'LineStyle', '--', 'LineWidth',4)
        hold on

        plot(1:H+1, Pctl_IRF1_q10(1)*ones(size(1:H+1)), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2)
        set(gca,'FontSize',50)
        xlim([1 (H+1)])

        sName = ['Fig_AggSh',num2str(sh_id),'_Pctl',num2str(ii),'.pdf'];
        saveas(figure(n_agg+ii), [figDir '\' sName] );
        close all

    end

end

