%--------------------------------------------------------------------------
% Script to generate Figure A2 (densities: transformed vs. raw) 
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
xmin = 0;
xmax = 3;
xn = 301;
xgrid = linspace(xmin, xmax, xn);

nfVARSpec = '10tc';
periods = unique(earnings_t);
T = length(periods);
K = 10;

sNameDir = ['fVAR',nfVARSpec];
loadDir = [pwd, '/', 'results' ,'/', sNameDir,'/'];
saveDir = [pwd, '/', 'figures' ,'/', sNameDir,'/'];
sNameFile = ['K',num2str(K),'_fVAR',nfVARSpec];    

% given period tt, draw the estimated density
PhatDensValue = csvread( [loadDir, sNameFile, '_PhatDensValue.csv'], 1, 0); % read table starting one row below
   
figure(1);clf;
set(figure(1),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(xgrid,PhatDensValue(1,:),'Color','b','LineStyle','-','LineWidth',2)
hold on

for tt = 2:T
    plot(xgrid,PhatDensValue(tt,:),'Color','b','LineStyle','-','LineWidth',2)
    hold on
end
       
set(gca,'FontSize',30)
xlim([0 4])
ylim([0 1])

sNameFig = '_Figure_A2_transf.pdf';
saveas(figure(1), [saveDir  sNameFile sNameFig] );
% close all
    
theta_sinh = 1.0;
ygrid = 1/(2*theta_sinh)*(exp(theta_sinh*xgrid) - exp(-theta_sinh*xgrid));
Jacobian = 1/2*(exp(theta_sinh*xgrid) + exp(-theta_sinh*xgrid));

figure(2);clf;
set(figure(2),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(ygrid,PhatDensValue(1,:)./Jacobian,'Color','b','LineStyle','-','LineWidth',2)
hold on    

for tt = 2:T
    plot(ygrid,PhatDensValue(tt,:)./Jacobian,'Color','b','LineStyle','-','LineWidth',2)
    hold on
end
 
set(gca,'FontSize',30)
xlim([0 4])
ylim([0 1])

sNameFig = '_Figure_A2_raw.pdf';
saveas(figure(2), [saveDir  sNameFile sNameFig] );
% close all

xgrid = linspace(0.01, 10, 100);
logx = log(xgrid);
ttheta = 1;
invsinhx = log(ttheta*xgrid + (ttheta^2*(xgrid.^2) + 1).^(1/2))/ttheta;


figure(3);clf;
set(figure(3),'PaperType','usletter','PaperOrientation','Landscape','PaperPosition',[0.1 0.1 11 8.5]);
plot(xgrid,invsinhx,'Color','b','LineStyle','-','LineWidth',3)
hold on    
plot(xgrid,logx,'Color','r','LineStyle','--','LineWidth',3)
hold on
plot(xgrid,xgrid,'Color','m','LineStyle',':','LineWidth',3)

set(gca,'FontSize',30)
xlim([0 4])

sNameFig = '_Figure_A2_invsinhy.pdf';
saveas(figure(3), [saveDir  sNameFile sNameFig] );


