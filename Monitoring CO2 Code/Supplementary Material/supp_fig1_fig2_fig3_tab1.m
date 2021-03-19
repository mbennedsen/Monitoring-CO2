%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "Designing a statistical procedure for monitoring global 
% carbon dioxide emissions" (2021) by Mikkel Bennedsen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code will produce Figures 1-3 and Table 1 from the Supplementary Material. Before running the code, please download
% the relevant data (freely) from https://doi.org/10.18160/gcp-2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Mikkel Bennedsen, February 2021.
% Code can be freely used and distributed. Please cite Bennedsen (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; close all;
%% init
str = {'Global_Carbon_Budget_2017v1.3.xlsx','Global_Carbon_Budget_2018v1.0.xlsx','Global_Carbon_Budget_2019v1.0.xlsx','Global_Carbon_Budget_2020v1.0.xlsx'};
pl_str = {'GCB2017','GCB2018','GCB2019','GCB2020'};

conc_1750 = 278;%276.5; %ppm
conc_1959 = 315.39; %ppm

maxp = 5;
maxq = 5;

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
    
    fig90 = figure(90);
    subplot(length(str),2,2*(iS-1)+1);
    autocorr(B_IM);
    title(['ACF: ',pl_str{iS}]);
    ylabel('Sample ACF');
    
    subplot(length(str),2,2*iS);
    parcorr(B_IM);
    title(['PACF: ',pl_str{iS}]);
    ylabel('Sample PACF');
    
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
    
    
    %% Disp
    [h_kpss,pValue_kpss] = kpsstest(B_IM,'lags',0:5, 'trend',false);%kpsstest(B_IM);
    [h_adf,pValue_adf] = adftest(B_IM,'lags',0:5);
    disp(['KPSS: h = ',num2str(h_kpss)]);%,', p-value = ',num2str(pValue)]);
    disp(['p-value = ',num2str(pValue_kpss)]);
    
    disp(['ADF: h = ',num2str(h_adf)]);%,', p-value = ',num2str(pValue)]);
    disp(['p-value = ',num2str(pValue_adf)]);
    
    [h_ks,p_ks,ksstat,cv_ks] = kstest(B_IM);
    [h_ad,p_ad,adstat,cv_ad] = adtest(B_IM);
    %% Fit AR(1) to the process; calc residuals
    y_fitted = c_ar1*B_IM(1:end-1);
    AR1_res = y_fitted - B_IM(2:end);
    
    sig2_hat = AR1_res'*AR1_res/(length(AR1_res)-1);

    allRES(1:length(AR1_res),iS) = AR1_res;
    
    output_BIM = [output_BIM;[length(B_IM),mean(B_IM),std(B_IM),skewness(B_IM),kurtosis(B_IM),c_ar1,sqrt(sig2_hat),NN,ksstat,adstat,DW,LB1,LB5]];
    
    
    fig91 = figure(91);
    subplot(length(str),2,2*(iS-1)+1);
    autocorr(AR1_res);
    title(['ACF, AR1 residuals: ',pl_str{iS}]);
    ylabel('Sample ACF');
    
    subplot(length(str),2,2*iS);
    parcorr(AR1_res);
    title(['PACF, AR1 residuals: ',pl_str{iS}]);
    ylabel('Sample PACF');
    
    %% QQ plot
    fig5 = figure(5);
    subplot(2,2,iS);
    qqplot(AR1_res);
    title(pl_str{iS});
    
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
    %% Disp
    [h_kpss,pValue_kpss] = kpsstest(AR1_res,'lags',0:5, 'trend',false);
    [h_adf,pValue_adf] = adftest(AR1_res,'lags',0:5);
    disp(['KPSS: h = ',num2str(h_kpss)]);%,', p-value = ',num2str(pValue)]);
    disp(['p-value = ',num2str(pValue_kpss)]);
    
    disp(['ADF: h = ',num2str(h_adf)]);%,', p-value = ',num2str(pValue)]);
    disp(['p-value = ',num2str(pValue_adf)]);

    %% Find best ARMA(p,q) model using BIC 
    LOGL = zeros(maxp+1,maxq+1); % Initialize
    PQ = zeros(maxp+1,maxq+1);
    for p = 0:maxp
        for q = 0:maxq
            mod = arima(p,0,q);
            [fit,~,logL] = estimate(mod,B_IM,'Display','off');
            LOGL(p+1,q+1) = logL;
            PQ(p+1,q+1) = p+q;
         end
    end

    LOGL = reshape(LOGL,(maxp+1)*(maxq+1),1);
    PQ = reshape(PQ,(maxp+1)*(maxq+1),1);
    [~,bic] = aicbic(LOGL,PQ+1,length(B_IM));

    BIC_output = reshape(bic,maxp+1,maxq+1);


    [val,indx] = min(BIC_output);
    [val2,indx2] = min(val);
    bic_p = indx(indx2)-1;
    bic_q = indx2-1;

    disp(' ');
    disp(['Optimal BIC: [p,q] = ',num2str([bic_p,bic_q])]);
end





