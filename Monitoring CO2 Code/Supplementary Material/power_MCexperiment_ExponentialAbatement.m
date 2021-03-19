%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "Designing a statistical procedure for monitoring global 
% carbon dioxide emissions" (2021) by Mikkel Bennedsen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Exponential Abatement: Monte Carlo power experiments.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Mikkel Bennedsen, February 2021.
% Code can be freely used and distributed. Please cite Bennedsen (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
rng(666);
%% Init
K = 60;  % Initial period
T = 30; % Number of periods to monitor (all periods = T+K)

DGP = 4;
rho = 0.90;
breaktime = 1;
m_vec = 0.05:0.05:1;

E_ANT = 9.9456;
gg = -0.0692;

MC = 10000;

if DGP == 1
    phi = 0.35; % AR coeff for BIM
    sig = 0.72; % Std. dev. for BIM innovations
    
    randomize_coeff = 0; % if 1, randomize phi; if 2 randomize both phi and sig at each revision.
elseif DGP == 2
    phi = 0.35/2; % AR coeff for BIM
    sig = 0.72; % Std. dev. for BIM innovations
    
    randomize_coeff = 0; % if 1, randomize phi; if 2 randomize both phi and sig at each revision.    
elseif DGP == 3
    phi = 0.35; % AR coeff for BIM
    sig = 0.72/2; % Std. dev. for BIM innovations
    
    randomize_coeff = 0; % if 1, randomize phi; if 2 randomize both phi and sig at each revision. 
elseif DGP == 4    
    lambda = 0.72;
    randomize_coeff = 2; % if 1, randomize phi; if 2 randomize both phi and sig at each revision.
else
    error('DGP not implemented');
end

%% Provide boundary function
g_fct = @(t,c)( c*sqrt(t) );
load(['crit_val_T',num2str(T),'_v01.mat']);

%% Calc critical boundaries
boundary_cond  = nan(T,3);
boundary_cond2 = nan(T,3);
Z_hat = nan(T,2);
for i = 1:T
    r = i;
    
    %%% Two-sided
    cst1 =  c_tilde_twoSide(1); % 5%
    boundary_cond(i,1) =  g_fct(r,cst1);
    
    cst2 =  c_tilde_twoSide(2); %10%
    boundary_cond(i,2) =  g_fct(r,cst2);
    
    cst3 = c_tilde_twoSide(3); % 32%
    boundary_cond(i,3) =  g_fct(r,cst3);
    
    %%% One-sided
    cst1 =  c_tilde_oneSide(1); % 5%
    boundary_cond2(i,1) =  -1*g_fct(r,cst1);
    
    cst2 =  c_tilde_oneSide(2); %10%
    boundary_cond2(i,2) =  -1*g_fct(r,cst2);
    
    cst3 = c_tilde_oneSide(3); % 32%
    boundary_cond2(i,3) =  -1*g_fct(r,cst3);
end
%% Run MC experiment
twoSide_res  = [];
oneSide_res  = [];
twoSide_det  = [];
oneSide_det  = [];
for iM = 1:length(m_vec)
    disp(iM/length(m_vec));
    
    %m = gg*(1-m_vec(iM));
    f = m_vec(iM);
    
    tau = K + breaktime;
    
    %% Simulation study
    detect_year = nan(MC,3);
    detect_year2 = nan(MC,3);
    detect_dummy = zeros(MC,3);
    detect_dummy2 = zeros(MC,3);
    for iMC = 1:MC
        t = 0; % Current period
        detect_period = -1;
        
        phi_save = 0;
        sig_save = 0;
        LRV_save = 0;
        
        while detect_period < 0 && t < T
            t = t + 1;
            
            
            %% Draw parameters
            if randomize_coeff == 1
                phi = unifrnd(0,1);
            elseif randomize_coeff == 2
                phi = unifrnd(0,1);
                sig = exprnd(lambda);
            end
            
            %% Simulate AR(1) process
            X = nan(K+t,1);
            X_H1 = nan(K+t,1);
            x0 = sig/sqrt(1-phi^2)*randn;
            if t == 1
                eps_t = randn(K+t,1);
            else
                eps_t = [rho*eps_t + sqrt(1-rho^2)*randn(K+t-1,1);
                    randn];
            end
            
            t_break = 0;
            for it = 1:(K+t)
                if it == 1
                    X(1) = phi*x0 + sig*eps_t(1);
                else
                    X(it) = phi*X(it-1) + sig*eps_t(it);
                end
                
                if it>tau-1
                    t_break = t_break+1;
                    E_path_act = (1-f)*(1+gg)^t_break*E_ANT(end) + f*E_ANT(end);
                    E_path_ast = (1+gg)^t_break*E_ANT(end);
                    
                    X_H1(it) = X(it) + E_path_ast - E_path_act;
                else
                    X_H1(it) = X(it);
                end
            end
            
            X = X_H1;
            %% Estimate AR params
            xxx = X(1:K-1); yyy = X(2:K);
            phi_hat = (xxx'*xxx)\xxx'*yyy;
            sig2_hat = sum( (yyy - xxx*phi_hat).^2)/(K-1);
            omega2_hat = sig2_hat/(1-phi_hat)^2;
            
            %% Estimate error process
            eps_hat = (X(K+t) - phi_hat*X(K+t-1))/sqrt(sig2_hat);
            
            %% Calculate detection stat
            if t == 1
                Ztilde = eps_hat;
            else
                Ztilde = Ztilde + eps_hat;
            end
            
            for iC = 1:3
                if Ztilde > boundary_cond(t,iC) || Ztilde < -boundary_cond(t,iC)
                    if iC == 1
                        detect_period = t;
                    end
                    if isnan(detect_year(iMC,iC)) == 1
                        detect_year(iMC,iC) = t;
                        detect_dummy(iMC,iC) = 1;
                    end
                    
                end
                
                if  Ztilde < boundary_cond2(t,iC)
                    if isnan(detect_year2(iMC,iC)) == 1
                        detect_year2(iMC,iC) = t;
                        detect_dummy2(iMC,iC) = 1;
                    end
                end
            end
            
        end
   
    end
    
    twoSide_res = [twoSide_res;nanmean(detect_dummy)];
    oneSide_res = [oneSide_res;nanmean(detect_dummy2)];
    
    twoSide_det = [twoSide_det;nanmean(detect_year)];
    oneSide_det = [oneSide_det;nanmean(detect_year2)];

end

%% Save
save(['Files/power_expAbate_DGP',num2str(DGP),'_K',num2str(K),'_T',num2str(T),'_v01'],'twoSide_res','oneSide_res','oneSide_det','MC');

