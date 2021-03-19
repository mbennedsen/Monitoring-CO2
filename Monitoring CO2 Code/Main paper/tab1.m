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

allBIM = nan(61,length(str));
allRES = nan(60,length(str));
output_BIM = [];
output_res = [];
output_rev = [];
for iS = 1:length(str)

    %% Load data
    dat = xlsread(str{iS},2);

    %% Create variables (notation as in Le Quere et al., 2017, see, e.g., Equation (1) or Table 2).
    if iS == 1
        disp('Analysing GCB2017 data...')
        t       = dat(:,1);
        B_IM    = dat(:,7);
        
        allBIM(1:length(B_IM),iS) = B_IM;
    elseif iS == 2        
        disp('Analysing GCB2018 data...')
        t       = dat(:,1);
        B_IM    = dat(:,7);
        
        allBIM(1:length(B_IM),iS) = B_IM;
        
    elseif iS == 3
        disp('Analysing GCB2019 data...')
        t       = dat(:,1+1);
        B_IM    = dat(:,7+1);
        
        allBIM(1:length(B_IM),iS) = B_IM;
    elseif iS == 4   
        disp('Analysing GCB2020 data...')
        t       = dat(:,1);
        B_IM    = dat(:,8);

        allBIM(1:length(B_IM),iS) = B_IM;
    else
        asff;
    end
    
    %% Make stuff
    v_t = B_IM;

    %% Calc N, DW
    n = length(v_t);

    m    = nan(4,1);
    m(1) = mean(v_t);
    for i = 2:4
        m(i) = mean( (v_t-m(1)).^i );
    end

    S = m(3)/sqrt(m(2)^3);
    K = m(4)/m(2)^2;
    NN = n*(S^2/6 + (K-3)^2/24);

    DW = sum( diff(v_t).^2 )/sum(v_t.^2);

    
    c_ar1 = mvregress(B_IM(1:end-1),B_IM(2:end));
    [h1,pValue1,LB1,cValue1] = lbqtest(B_IM,'Lags',1);
    [h2,pValue2,LB5,cValue2] = lbqtest(B_IM,'Lags',5);
    
   
    
    [h_ks,p_ks,ksstat,cv_ks] = kstest(B_IM);
    [h_ad,p_ad,adstat,cv_ad] = adtest(B_IM);
    %% Fit AR(1) to the process; calc residuals
    y_fitted = c_ar1*B_IM(1:end-1);
    AR1_res = y_fitted - B_IM(2:end);
    
    sig2_hat = AR1_res'*AR1_res/(length(AR1_res)-1);

    allRES(1:length(AR1_res),iS) = AR1_res;
    
    output_BIM = [output_BIM;[length(B_IM),mean(B_IM),std(B_IM),skewness(B_IM),kurtosis(B_IM),c_ar1,sqrt(sig2_hat),NN,ksstat,adstat,DW,LB1,LB5]];
    
    %% Make stuff       
    AR1_res = AR1_res/sqrt(sig2_hat);
    
    v_t = AR1_res;
    %% Calc N, DW for fitted AR(1) residuals
    n = length(v_t);

    m    = nan(4,1);
    m(1) = mean(v_t);
    for i = 2:4
        m(i) = mean( (v_t-m(1)).^i );
    end

    S = m(3)/sqrt(m(2)^3);
    K = m(4)/m(2)^2;
    NN = n*(S^2/6 + (K-3)^2/24);

    DW = sum( diff(v_t).^2 )/sum(v_t.^2);
    
    c_ar1_0 = mvregress(AR1_res(1:end-1),AR1_res(2:end));
    [h1,pValue1,LB1,cValue1] = lbqtest(AR1_res,'Lags',1);
    [h2,pValue2,LB5,cValue2] = lbqtest(AR1_res,'Lags',5);
    
    y_fitted2 = c_ar1_0*AR1_res(1:end-1);
    AR1_res2 = y_fitted2 - AR1_res(2:end);
    
    sig2_hat2 = AR1_res2'*AR1_res2/(length(AR1_res2)-1);
    
    
    [h_ks,p_ks,ksstat,cv_ks] = kstest(AR1_res);
    [h_ad,p_ad,adstat,cv_ad] = adtest(AR1_res);
    
    output_res = [output_res;[length(AR1_res),mean(AR1_res),std(AR1_res),skewness(AR1_res),kurtosis(AR1_res),c_ar1_0,sqrt(sig2_hat2),NN,ksstat,adstat,DW,LB1,LB5]];
 
end

%% Print to screen

disp(' ');
disp('----------------------- Diagnostics - BIM ----------------------- ');
disp('    numObs     Mean      Std.      Skew.     Kurt.     AR(1)      sig       JB        KS        AD        DW       LB(1)     LB(5)');
disp(round(output_BIM,2));

disp(' ');
disp('----------------------- Diagnostics - Residuals ----------------------- ');
disp('    numObs     Mean      Std.      Skew.     Kurt.     AR(1)      sig       JB        KS        AD        DW       LB(1)     LB(5)');
disp(round(output_res,2));




