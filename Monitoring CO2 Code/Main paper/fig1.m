%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "Designing a statistical procedure for monitoring global 
% carbon dioxide emissions" (2021) by Mikkel Bennedsen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code will produce Figure 1 from the paper. Before running the code, please download
% the relevant data (freely) from https://doi.org/10.18160/gcp-2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Mikkel Bennedsen, February 2021.
% Code can be freely used and distributed. Please cite Bennedsen (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;

%% init

str = {'Global_Carbon_Budget_2017v1.3.xlsx','Global_Carbon_Budget_2018v1.0.xlsx','Global_Carbon_Budget_2019v1.0.xlsx','Global_Carbon_Budget_2020v1.0.xlsx'};

conc_1750 = 278;%276.5; %ppm
conc_1959 = 315.39; %ppm

maxp = 5;
maxq = 5;

allBIM = nan(60,3);
output_BIM = [];
output_res = [];
output_rev = [];

fig = figure(1);
for iS = 1:length(str)

    %% Load data
    dat = xlsread(str{iS},2);

    %% Create variables (notation as in Le Quere et al., 2017, see, e.g., Equation (1) or Table 2).
    if iS == 1
        disp('Analysing GCB2017 data...')
        t       = dat(:,1);
        E_FF    = dat(:,2);
        E_LUC   = dat(:,3);
        G_ATM   = dat(:,4);
        S_OCEAN = dat(:,5);
        S_LAND  = dat(:,6);
        B_IM    = dat(:,7);

        E_ANT = E_FF + E_LUC;


    elseif iS == 2

        disp('Analysing GCB2018 data...')
        t       = dat(:,1);
        E_FF    = dat(:,2);
        E_LUC   = dat(:,3);
        G_ATM   = dat(:,4);
        S_OCEAN = dat(:,5);
        S_LAND  = dat(:,6);
        B_IM    = dat(:,7);

        E_ANT = E_FF + E_LUC;

        
        
    elseif iS == 3
        disp('Analysing GCB2019 data...')
        t       = dat(:,1+1);
        E_FF    = dat(:,2+1);
        E_LUC   = dat(:,3+1);
        G_ATM   = dat(:,4+1);
        S_OCEAN = dat(:,5+1);
        S_LAND  = dat(:,6+1);
        B_IM    = dat(:,7+1);

        E_ANT = E_FF + E_LUC;
        
     elseif iS == 4
        disp('Analysing GCB2020 data...')
        t        = dat(:,1);
        E_FF     = dat(:,2);
        E_LUC    = dat(:,3);
        G_ATM    = dat(:,4);
        S_OCEAN  = dat(:,5);
        S_LAND   = dat(:,6);
        S_CEMENT = dat(:,7);
        B_IM     = dat(:,8);

        E_ANT = E_FF + E_LUC;
        
        E_FF = E_FF-S_CEMENT; % Include cement carbonation sink into E_FF

    else
        asff;
    end
    
    %% plot
    
    
    subplot(3,2,1);
    if iS == 1
        cols = [0, 0.4470, 0.7410];
    elseif iS == 2
        cols = [0.8500, 0.3250, 0.0980];
    elseif iS == 3
        cols = [0.9290, 0.6940, 0.1250];
    elseif iS == 4
        cols = [0.4940, 0.1840, 0.5560];
    else
        asdf;
    end
        
    plot(t,E_FF,'LineWidth',1.5,'Color',cols), hold on
    plot(t,E_LUC,'LineWidth',1.5,'Color',cols), hold on
    ylabel('GtC/yr','Interpreter','latex','FontSize',10);
    title('Anthropogenic emissions','Interpreter','latex','FontSize',10);
    axis([1959,2020,0,10.5]);

    subplot(3,2,2);
    plot(t,G_ATM,'LineWidth',1.5), hold on
    ylabel('GtC/yr','Interpreter','latex','FontSize',10);
    title('Atmospheric growth','Interpreter','latex','FontSize',10);

    subplot(3,2,3);
    plot(t,S_OCEAN,'LineWidth',1.5), hold on
    ylabel('GtC/yr','Interpreter','latex','FontSize',10);
    title('Ocean sink','Interpreter','latex','FontSize',10);
    if iS == 4
        leg1 = legend('GCB2017','GCB2018','GCB2019','GCB2020','Location','Best');
        set(leg1,'Interpreter','latex');
        set(leg1,'FontSize',8);
        
        legend('boxoff');
    end

    subplot(3,2,4);
    plot(t,S_LAND,'LineWidth',1.5), hold on
    ylabel('GtC/yr','Interpreter','latex','FontSize',10);
    title('Land sink','Interpreter','latex','FontSize',10);

    subplot(3,2,[5,6]);
    if iS == 1
        plot([t;t(end)+(1:5)'],zeros(length(t)+5,1),'k-'), hold on
    end
    plot(t,B_IM,'LineWidth',1.5), hold on
    ylabel('GtC/yr','Interpreter','latex','FontSize',10);

    axis([1959,2020,-2.45,2.45]);
    title('Budget imbalance','Interpreter','latex','FontSize',10);


end




