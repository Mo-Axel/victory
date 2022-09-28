%--------------------------------------------------------------------------
% Script to generate Figure on MDD 
%--------------------------------------------------------------------------

clear; 
clc;
close all

set(0,'defaultTextInterpreter','latex');

% set specs 
xmin = 0;
xmax = 8;
xn = 100;
xgrid = linspace(xmin, xmax, xn);

% colsel does not count first three columns (hyperparameters)
colsel = 7; 
K = 20;

nfVARSpec = '10tc';
%7tc;
nModSpec  = '1';
nMCMCSpec = '1';
sName     = ['fVAR', nfVARSpec, '_MDD', nModSpec, '_MCMC', nMCMCSpec];

figsaveDir= [pwd, '\', 'Figures' ,'\', sName,'\'];
[~, ~, ~]  = mkdir(figsaveDir);

% load MDDs: use Bayes vs Laplace
mddDir = [pwd, '\', 'Results' ,'\', sName, '\'];
lambda_MDD = readmatrix([mddDir, sName, '_MDD_Laplace_sum.csv']);
%, 1, 0
lambda1 = lambda_MDD(:,1);
lambda2 = lambda_MDD(:,2);
lambda3 = lambda_MDD(:,3);
MDD = lambda_MDD(:,3+colsel);


maxMDD_id = find(MDD == max(MDD));

maxlambda1 = lambda1(maxMDD_id);
maxlambda1_id = find(lambda1 == maxlambda1);

maxlambda2 = lambda2(maxMDD_id);
maxlambda3 = lambda3(maxMDD_id);

% fix lambda1 = maxlambda1, get the surface plot for MDD

grid_lambda1 = unique(lambda1);
grid_lambda2 = unique(lambda2);
grid_lambda3 = unique(lambda3);

Ngrid_lambda1 = length(grid_lambda1);
Ngrid_lambda2 = length(grid_lambda2);
Ngrid_lambda3 = length(grid_lambda3);

[X, Y]=meshgrid(unique(lambda3),unique(lambda2)); % be careful with the order

for ii = maxlambda1_id(1):maxlambda1_id(end)
 surfaceMDD(find(grid_lambda2 == lambda2(ii)),find(grid_lambda3 == lambda3(ii))) = MDD(ii);
end

figure(1);clf;
set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
surf(X,Y,surfaceMDD)
set(gca,'FontSize',30)
%xlabel('\lambda_3')
%ylabel('\lambda_2')
sNameFile = ['K', num2str(K),'_Figure_MDD3d.pdf'];    
saveas(figure(1), [figsaveDir sNameFile] );
close all


% at maximum MDD's lambda1 and labmda2 value
maxlambda12_id = find(lambda1 == maxlambda1 & lambda2 == maxlambda2);
MDD_maxlambda12 = MDD(maxlambda12_id);
MDD_max_id = find(MDD_maxlambda12 == max(MDD_maxlambda12));
MDD_norm = MDD_maxlambda12(MDD_max_id);

figure(2);clf; 
set(figure(2),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(log(grid_lambda3), MDD_maxlambda12 - MDD_norm,'Color','b','LineStyle','-','LineWidth',4)
xlim([log(grid_lambda3(1)) log(grid_lambda3(end))])
%plot(1:Ngrid_lambda3, MDD(maxlambda12_id),'Color','b','LineStyle','-','LineWidth',4)
%xtick_locs = [1 round(Ngrid_lambda3/3) round(Ngrid_lambda3*2/3) Ngrid_lambda3];
%xticks(xtick_locs); 
%xticklabels({num2str(grid_lambda3(xtick_locs(1))), num2str(grid_lambda3(xtick_locs(2))), ...
%    num2str(grid_lambda3(xtick_locs(3))), num2str(grid_lambda3(xtick_locs(4)))})
set(gca,'FontSize',50)

sNameFile = ['K', num2str(K),'_Figure_MDD2d_lambda3.pdf'];    
saveas(figure(2), [figsaveDir sNameFile] );
close all


% at maximum MDD's lambda1 and labmda3 value
maxlambda13_id = find(lambda1 == maxlambda1 & lambda3 == maxlambda3);
MDD_maxlambda13 = MDD(maxlambda13_id);
MDD_max_id = find(MDD_maxlambda13 == max(MDD_maxlambda13));
MDD_norm = MDD_maxlambda13(MDD_max_id);

figure(3);clf; 
set(figure(3),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(log(grid_lambda2), MDD_maxlambda13 - MDD_norm,'Color','b','LineStyle','-','LineWidth',4)
xlim([log(grid_lambda2(1)) log(grid_lambda2(end))])
%plot(1:Ngrid_lambda2, MDD(maxlambda13_id),'Color','b','LineStyle','-','LineWidth',4)
%xtick_locs = [1 round(Ngrid_lambda2/3) round(Ngrid_lambda2*2/3) Ngrid_lambda2];
%xticks(xtick_locs); 
%xticklabels({num2str(grid_lambda2(xtick_locs(1))), num2str(grid_lambda2(xtick_locs(2))), ...
%    num2str(grid_lambda2(xtick_locs(3))), num2str(grid_lambda2(xtick_locs(4)))})
set(gca,'FontSize',50)

sNameFile = ['K', num2str(K),'_Figure_MDD2d_lambda2.pdf'];    
saveas(figure(3), [figsaveDir sNameFile] );
close all

% at maximum MDD's lambda2 and labmda3 value
maxlambda23_id = find(lambda2 == maxlambda2 & lambda3 == maxlambda3);
MDD_maxlambda23 = MDD(maxlambda23_id);
MDD_max_id = find(MDD_maxlambda23 == max(MDD_maxlambda23));
MDD_norm = MDD_maxlambda23(MDD_max_id);

figure(4);clf; 
set(figure(4),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(log(grid_lambda1), MDD_maxlambda23 - MDD_norm, 'Color','b','LineStyle','-','LineWidth',4)
xlim([log(grid_lambda1(1)) log(grid_lambda1(end))])
%plot(1:Ngrid_lambda1, MDD(maxlambda23_id),'Color','b','LineStyle','-','LineWidth',4)
%xtick_locs = [1 round(Ngrid_lambda1/3) round(Ngrid_lambda1*2/3) Ngrid_lambda1];
%xticks(xtick_locs); 
%xticklabels({num2str(grid_lambda1(xtick_locs(1))), num2str(grid_lambda1(xtick_locs(2))), ...
%    num2str(grid_lambda1(xtick_locs(3))), num2str(grid_lambda1(xtick_locs(4)))})
set(gca,'FontSize',50)

sNameFile = ['K', num2str(K),'_Figure_MDD2d_lambda1.pdf'];    
saveas(figure(4), [figsaveDir sNameFile] );
close all


