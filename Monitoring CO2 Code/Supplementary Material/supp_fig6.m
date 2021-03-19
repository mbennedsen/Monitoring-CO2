%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "Designing a statistical procedure for monitoring global 
% carbon dioxide emissions" (2021) by Mikkel Bennedsen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code will produce Figure 6 from the Supplementary Material..
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Mikkel Bennedsen, February 2021.
% Code can be freely used and distributed. Please cite Bennedsen (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
%% Init
T = 30; % Number of periods to monitor
m_vec1 = [0.10,0.20,0.30]; % Exp
gg = -0.0692;
m_vec2 = [0.10,0.20,0.30];% Lin

E_ANT = 9.9456;
E_ast_slope = -(0.5*9.0430-E_ANT)/11;
%% Calc
E_paths_exp = nan(T+1,length(m_vec1),2);
for iM = 1:length(m_vec1)
    %m = gg*(1-m_vec1(iM));
    f = m_vec1(iM);
    for iT = 0:T
        
        E_path_act = (1-f)*(1+gg)^iT*E_ANT(end) + f*E_ANT(end);
        E_path_ast = (1+gg)^iT*E_ANT(end);
        
        E_paths_exp(iT+1,iM,1) = E_path_act;
        E_paths_exp(iT+1,iM,2) = E_path_ast;
        
    end
    
end

E_paths_lin = nan(T+1,length(m_vec2),2);
for iM = 1:length(m_vec2)
    f = m_vec2(iM);

    for iT = 0:T
        
        E_path_act = (1-f)*max(E_ANT(end) - E_ast_slope*iT,0) + f*E_ANT(end);
        E_path_ast = E_ANT(end) - E_ast_slope*iT;
        
        E_paths_lin(iT+1,iM,1) = max(E_path_act,0);
        E_paths_lin(iT+1,iM,2) = max(E_path_ast,0);
        
    end
    
end


%% plot
fig1 = figure;
subplot(3,2,2)
plot(0:T,E_paths_lin(:,1,1),'b-','LineWidth',1.5), hold on
plot(0:T,E_paths_lin(:,1,2),'r-.','LineWidth',1.5), hold on
ylabel('GtC/yr','Interpreter','latex','FontSize',10);
title(['Linear abatement ($m = ',num2str(m_vec2(1)),'$)'],'Interpreter','latex','FontSize',10);

subplot(3,2,4)
plot(0:T,E_paths_lin(:,2,1),'b-','LineWidth',1.5), hold on
plot(0:T,E_paths_lin(:,2,2),'r-.','LineWidth',1.5), hold on
ylabel('GtC/yr','Interpreter','latex','FontSize',10);
title(['Linear abatement ($m = ',num2str(m_vec2(2)),'$)'],'Interpreter','latex','FontSize',10);

subplot(3,2,6)
plot(0:T,E_paths_lin(:,3,1),'b-','LineWidth',1.5), hold on
plot(0:T,E_paths_lin(:,3,2),'r-.','LineWidth',1.5), hold on
xlabel('$t$','Interpreter','latex','FontSize',10);
ylabel('GtC/yr','Interpreter','latex','FontSize',10);
title(['Linear abatement ($m = ',num2str(m_vec2(3)),'$)'],'Interpreter','latex','FontSize',10);




subplot(3,2,1)
plot(0:T,E_paths_exp(:,1,1),'b-','LineWidth',1.5), hold on
plot(0:T,E_paths_exp(:,1,2),'r-.','LineWidth',1.5), hold on
ylabel('GtC/yr','Interpreter','latex','FontSize',10);
title(['Exponential abatement ($m = ',num2str(m_vec1(1)),'$)'],'Interpreter','latex','FontSize',10);

leg1 = legend('Actual emissions ($E_t^{FF}$)','Reported emissions ($E_t^{FF,*}$)','Location','NorthEast');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',8);
legend('boxoff');


subplot(3,2,3)
plot(0:T,E_paths_exp(:,2,1),'b-','LineWidth',1.5), hold on
plot(0:T,E_paths_exp(:,2,2),'r-.','LineWidth',1.5), hold on

ylabel('GtC/yr','Interpreter','latex','FontSize',10);
title(['Exponential abatement ($m = ',num2str(m_vec1(2)),'$)'],'Interpreter','latex','FontSize',10);

subplot(3,2,5)
plot(0:T,E_paths_exp(:,3,1),'b-','LineWidth',1.5), hold on
plot(0:T,E_paths_exp(:,3,2),'r-.','LineWidth',1.5), hold on
ylabel('GtC/yr','Interpreter','latex','FontSize',10);
xlabel('$t$','Interpreter','latex','FontSize',10);
title(['Exponential abatement ($m = ',num2str(m_vec1(3)),'$)'],'Interpreter','latex','FontSize',10);
