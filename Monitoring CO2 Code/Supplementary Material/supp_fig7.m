%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "Designing a statistical procedure for monitoring global 
% carbon dioxide emissions" (2021) by Mikkel Bennedsen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code will produce Figure 7 from the Supplementary Material. Power experiments for the Exponential Abatement 
% case can be conducted by running the file "power_MCexperiment_ExponentialAbatement.m".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Mikkel Bennedsen, February 2021.
% Code can be freely used and distributed. Please cite Bennedsen (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
%% init
DGP1 = 1;
DGP2 = 2;
DGP3 = 3;
DGP4 = 4;

K = 60;
T = 30;
m_vec = 0.05:0.05:1;

%% Load
load(['Files/power_expAbate_DGP',num2str(DGP1),'_K',num2str(K),'_T',num2str(T)]);

%% plot
fig1  = figure;
subplot(4,2,1);
plot(m_vec,oneSide_res(:,1),'-','LineWidth',1.5), hold on
plot(m_vec,oneSide_res(:,2),'--','LineWidth',1.5), hold on
plot(m_vec,oneSide_res(:,3),'-.','LineWidth',1.5), hold on
axis([0/100,m_vec(end),0,1]);
ylabel('Power','Interpreter','latex','FontSize',10);
title('DGP1','Interpreter','latex','FontSize',10);
leg1 = legend('$\alpha = 5\%$','$\alpha = 10\%$','$\alpha = 32\%$','Location','Best');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',8);
legend('boxoff');
%% plot
subplot(4,2,2);
plot(m_vec(1:end),oneSide_det(1:end,1),'-','LineWidth',1.5), hold on
plot(m_vec(1:end),oneSide_det(1:end,2),'--','LineWidth',1.5), hold on
plot(m_vec(1:end),oneSide_det(1:end,3),'-.','LineWidth',1.5), hold on
plot([0,m_vec],5*ones(length(m_vec)+1,1),'k--','LineWidth',1.0), hold on
axis([0/100,m_vec(end),0,17.5]);
ylabel('Average detection time','Interpreter','latex','FontSize',10);
title('DGP1','Interpreter','latex','FontSize',10);

%% Load
clear oneSide_res;
load(['Files/power_expAbate_DGP',num2str(DGP2),'_K',num2str(K),'_T',num2str(T)]);

%% plot
subplot(4,2,3);
plot(m_vec,oneSide_res(:,1),'-','LineWidth',1.5), hold on
plot(m_vec,oneSide_res(:,2),'--','LineWidth',1.5), hold on
plot(m_vec,oneSide_res(:,3),'-.','LineWidth',1.5), hold on
axis([0/100,m_vec(end),0,1]);
ylabel('Power','Interpreter','latex','FontSize',10);
title('DGP2','Interpreter','latex','FontSize',10);

%% plot
subplot(4,2,4);
plot(m_vec(1:end),oneSide_det(1:end,1),'-','LineWidth',1.5), hold on
plot(m_vec(1:end),oneSide_det(1:end,2),'--','LineWidth',1.5), hold on
plot(m_vec(1:end),oneSide_det(1:end,3),'-.','LineWidth',1.5), hold on
plot([0,m_vec],5*ones(length(m_vec)+1,1),'k--','LineWidth',1.0), hold on
axis([0/100,m_vec(end),0,17.5]);
ylabel('Average detection time','Interpreter','latex','FontSize',10);
title('DGP2','Interpreter','latex','FontSize',10);

%% Load
clear oneSide_res;
load(['Files/power_expAbate_DGP',num2str(DGP3),'_K',num2str(K),'_T',num2str(T)]);

%% plot
subplot(4,2,5);
plot(m_vec,oneSide_res(:,1),'-','LineWidth',1.5), hold on
plot(m_vec,oneSide_res(:,2),'--','LineWidth',1.5), hold on
plot(m_vec,oneSide_res(:,3),'-.','LineWidth',1.5), hold on
axis([0/100,m_vec(end),0,1]);
ylabel('Power','Interpreter','latex','FontSize',10);
title('DGP3','Interpreter','latex','FontSize',10);

%% plot
subplot(4,2,6);
plot(m_vec(1:end),oneSide_det(1:end,1),'-','LineWidth',1.5), hold on
plot(m_vec(1:end),oneSide_det(1:end,2),'--','LineWidth',1.5), hold on
plot(m_vec(1:end),oneSide_det(1:end,3),'-.','LineWidth',1.5), hold on
plot([0,m_vec],5*ones(length(m_vec)+1,1),'k--','LineWidth',1.0), hold on
axis([0/100,m_vec(end),0,17.5]);
ylabel('Average detection time','Interpreter','latex','FontSize',10);
title('DGP3','Interpreter','latex','FontSize',10);

%% Load
clear oneSide_res;
load(['Files/power_expAbate_DGP',num2str(DGP4),'_K',num2str(K),'_T',num2str(T)]);

%% plot
subplot(4,2,7);
plot(m_vec,oneSide_res(:,1),'-','LineWidth',1.5), hold on
plot(m_vec,oneSide_res(:,2),'--','LineWidth',1.5), hold on
plot(m_vec,oneSide_res(:,3),'-.','LineWidth',1.5), hold on
axis([0/100,m_vec(end),0,1]);
xlabel('$m$','Interpreter','latex','FontSize',10);
ylabel('Power','Interpreter','latex','FontSize',10);
title('DGP4','Interpreter','latex','FontSize',10);

%% plot
subplot(4,2,8);
plot(m_vec(1:end),oneSide_det(1:end,1),'-','LineWidth',1.5), hold on
plot(m_vec(1:end),oneSide_det(1:end,2),'--','LineWidth',1.5), hold on
plot(m_vec(1:end),oneSide_det(1:end,3),'-.','LineWidth',1.5), hold on
plot([0,m_vec],5*ones(length(m_vec)+1,1),'k--','LineWidth',1.0), hold on
axis([0/100,m_vec(end),0,17.5]);

xlabel('$m$','Interpreter','latex','FontSize',10);
ylabel('Average detection time','Interpreter','latex','FontSize',10);
title('DGP4','Interpreter','latex','FontSize',10);