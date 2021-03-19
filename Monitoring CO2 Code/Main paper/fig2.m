%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "Designing a statistical procedure for monitoring global 
% carbon dioxide emissions" (2021) by Mikkel Bennedsen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code will produce Figure 2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Mikkel Bennedsen, February 2021.
% Code can be freely used and distributed. Please cite Bennedsen (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
%% Init
addpath('Files');
Tvec = 30; % Number of periods to monitor (T)    

%% Stuff
Tmax = max(Tvec);
g_fct = @(t,c)( c*sqrt(t) ); %g-fct

boundary_cond  = nan(Tmax,3,3);
boundary_cond2 = nan(Tmax,3,3);
for j = 1:length(Tvec)
    T = Tvec(j);
    load(['crit_val_T',num2str(T),'_v01.mat']);

    for i = 1:T
        r = i;

        %%% Two-sided
        cst1 =  c_tilde_twoSide(1); % 5%
        boundary_cond(i,1,j) =  g_fct(r,cst1);

        cst2 =  c_tilde_twoSide(2); %10%
        boundary_cond(i,2,j) =  g_fct(r,cst2);

        cst3 = c_tilde_twoSide(3); % 32%
        boundary_cond(i,3,j) =  g_fct(r,cst3);

        %%% One-sided
        cst1 =  c_tilde_oneSide(1); % 5%
        boundary_cond2(i,1,j) =  g_fct(r,cst1);

        cst2 =  c_tilde_oneSide(2); %10%
        boundary_cond2(i,2,j) =  g_fct(r,cst2);

        cst3 = c_tilde_oneSide(3); % 32%
        boundary_cond2(i,3,j) =  g_fct(r,cst3);
    end
end


%% plot
fig = figure;

plot(1:Tmax,boundary_cond2(:,1,1),'b-','LineWidth',1.5), hold on
plot(1:Tmax,boundary_cond2(:,2,1),'b--','LineWidth',1.5), hold on
plot(1:Tmax,boundary_cond2(:,3,1),'b-.','LineWidth',1.5), hold on

xlabel('t','Interpreter','latex','FontSize',10);
ylabel('Critical value, $C_t^{\alpha}$','Interpreter','latex','FontSize',10);
leg1 = legend('$\alpha = 5\%$','$\alpha = 10\%$','$\alpha = 32\%$','Location','Best');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',8);
legend('boxoff');

%% Disp (Table 2)
round(boundary_cond2(1:10,:,1),2)'
