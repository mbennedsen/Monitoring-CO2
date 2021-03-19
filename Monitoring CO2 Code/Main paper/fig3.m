%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "Designing a statistical procedure for monitoring global 
% carbon dioxide emissions" (2021) by Mikkel Bennedsen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code will produce Figure 3.
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

%%
fig = figure;
subplot(1,3,1);
plot(2019+(0:T),E_paths_exp(:,1,1),'b-','LineWidth',1.5), hold on
plot(2019+(0:T),E_paths_exp(:,1,2),'r-.','LineWidth',1.5), hold on
ylabel('GtC/yr','Interpreter','latex','FontSize',10);
xlabel('$t$','Interpreter','latex','FontSize',10);
axis([2019,2050,0,11]);
title('$m = 0.10$','Interpreter','latex','FontSize',10);

subplot(1,3,2);
plot(2019+(0:T),E_paths_exp(:,2,1),'b-','LineWidth',1.5), hold on
plot(2019+(0:T),E_paths_exp(:,1,2),'r-.','LineWidth',1.5), hold on
ylabel('GtC/yr','Interpreter','latex','FontSize',10);
xlabel('$t$','Interpreter','latex','FontSize',10);
axis([2019,2050,0,11]);
title('$m = 0.20$','Interpreter','latex','FontSize',10);

subplot(1,3,3);
plot(2019+(0:T),E_paths_exp(:,3,1),'b-','LineWidth',1.5), hold on
plot(2019+(0:T),E_paths_exp(:,1,2),'r-.','LineWidth',1.5), hold on
ylabel('GtC/yr','Interpreter','latex','FontSize',10);
xlabel('$t$','Interpreter','latex','FontSize',10);
axis([2019,2050,0,11]);
title('$m = 0.30$','Interpreter','latex','FontSize',10);

