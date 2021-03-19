%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "Designing a statistical procedure for monitoring global 
% carbon dioxide emissions" (2021) by Mikkel Bennedsen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code will produce Tables 2-3 from the Supplementary Material. Before running the code, please download
% the relevant data (freely) from https://doi.org/10.18160/gcp-2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Mikkel Bennedsen, February 2021.
% Code can be freely used and distributed. Please cite Bennedsen (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;

%% init
addpath('Files');

pl_str = {'GCB2017','GCB2018','GCB2019','GCB2020'};
str = {'Global_Carbon_Budget_2017v1.3.xlsx','Global_Carbon_Budget_2018v1.0.xlsx','Global_Carbon_Budget_2019v1.0.xlsx','Global_Carbon_Budget_2020v1.0.xlsx'};
pl_cl = {'r-','c-.','g--','b-'};
conc_1750 = 278;%276.5; %ppm
conc_1959 = 315.39; %ppm

mm_vals = [0,-0.05,-0.10]; % Magnitude of misreporting.
tau_vals = [0,0,0]; % Time to introduce misreporting.

cusum_p2 = [];
output = nan(length(str),5,length(mm_vals));
for iS = 1:length(str)

    %% Load data
    dat = xlsread(str{iS},2);
    
    for iV = 1:length(mm_vals)
        mm = mm_vals(iV);
        tau = tau_vals(iV);

        %% Create variables (notation as in Le Quere et al., 2017, see, e.g., Equation (1) or Table 2).
        if iS == 1
            disp('Analysing GCB2017 data...')
            t       = dat(:,1);
            E_FF    = dat(:,2);
            E_LUC   = dat(:,3);
            G_ATM   = dat(:,4);
            S_OCEAN = dat(:,5);
            S_LAND  = dat(:,6);
            B_IM_tr    = dat(:,7);

            mm_vec = [zeros(tau,1);
                      mm*ones(length(E_FF)-tau,1)];

            CB_IM = E_FF.*(1+mm_vec) + E_LUC - G_ATM - S_OCEAN - S_LAND;

        elseif iS == 2
            disp('Analysing GCB2018 data...')
            t       = dat(:,1);
            E_FF    = dat(:,2);
            E_LUC   = dat(:,3);
            G_ATM   = dat(:,4);
            S_OCEAN = dat(:,5);
            S_LAND  = dat(:,6);
            B_IM_tr    = dat(:,7);

            mm_vec = [zeros(tau,1);
                      mm*ones(length(E_FF)-tau,1)];

            CB_IM = E_FF.*(1+mm_vec) + E_LUC - G_ATM - S_OCEAN - S_LAND;

        elseif iS == 3
            disp('Analysing GCB2019 data...')
            t       = dat(:,1+1);
            E_FF    = dat(:,2+1);
            E_LUC   = dat(:,3+1);
            G_ATM   = dat(:,4+1);
            S_OCEAN = dat(:,5+1);
            S_LAND  = dat(:,6+1);
            B_IM_tr    = dat(:,7+1);

            mm_vec = [zeros(tau,1);
                      mm*ones(length(E_FF)-tau,1)];

            CB_IM = E_FF.*(1+mm_vec) + E_LUC - G_ATM - S_OCEAN - S_LAND;

        elseif iS == 4
            disp('Analysing GCB2020 data...')
            t       = dat(:,1);
            E_FF    = dat(:,2);
            E_LUC   = dat(:,3);
            G_ATM   = dat(:,4);
            S_OCEAN = dat(:,5);
            S_LAND  = dat(:,6);
            S_CEMENT = dat(:,7);
            B_IM_tr    = dat(:,8);

            mm_vec = [zeros(tau,1);
                      mm*ones(length(E_FF)-tau,1)];

            CB_IM = E_FF.*(1+mm_vec) + E_LUC - G_ATM - S_OCEAN - S_LAND - S_CEMENT;

            
        else
            asff;
        end




        %% Estimate parameters
        X = CB_IM;
        n = length(X);
        xxx = [ones(n-1,1),X(1:end-1)]; yyy = X(2:end); 
        beta_hat = (xxx'*xxx)\xxx'*yyy;
        phi_hat = beta_hat(2);
        eps_m_hat = yyy - xxx*beta_hat;
        sig2_hat = sum( (eps_m_hat).^2)/(n-1);

        omega2_hat = sig2_hat/(1-phi_hat)^2;
        
        %% Test for non-zero mean
        t_stat1 = mean(CB_IM)/(std(CB_IM)/sqrt(n));
        t_stat2 = mean(CB_IM)/sqrt(omega2_hat/n);
        
        output(iS,:,iV) = [n,mean(CB_IM),std(CB_IM),t_stat1,t_stat2];
    end
end

%% Output
disp(' ');
disp('Table 2 (Original data)');
disp(output(:,:,1));

disp(' ');
disp('Table 3 (m = 0.05)');
disp(output(:,:,2));

disp(' ');
disp('Table 3 (m = 0.10)');
disp(output(:,:,3));