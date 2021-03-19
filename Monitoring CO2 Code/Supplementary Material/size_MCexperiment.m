%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "Designing a statistical procedure for monitoring global 
% carbon dioxide emissions" (2021) by Mikkel Bennedsen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Monte Carlo size experiments.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Mikkel Bennedsen, February 2021.
% Code can be freely used and distributed. Please cite Bennedsen (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
rng(666);
%% Init
Kvec = [30,60,120];  % Initial period
Tvec = [30,60,120]; % Number of periods to monitor (all periods = T+K)

DGP = 1;
rho = 0.90;

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
    

%% Specify boundary function
g_fct = @(t,c)( c*sqrt(t) );
%% Get boundary functions
twoSide_res = [];
oneSide_res = [];
for iK = 1:length(Kvec)
    disp(iK/length(Kvec));
    
    K = Kvec(iK);
    
    for iT = 1:length(Tvec)
        T = Tvec(iT);
        
        load(['crit_val_T',num2str(T),'_v01.mat']);

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
        %% Simulation study
        detect_year = nan(MC,3);
        detect_year2 = nan(MC,3);
        detect_dummy = zeros(MC,3);
        detect_dummy2 = zeros(MC,3);
        for iMC = 1:MC
            x0 = 0;

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
                if t == 1
                    eps_t = randn(K+t,1);
                else
                    eps_t = [rho*eps_t + sqrt(1-rho^2)*randn(K+t-1,1);
                             randn];
                end
                for it = 1:(K+t)
                    if it == 1
                        X(1) = phi*x0 + sig*eps_t(1);
                    else
                        X(it) = phi*X(it-1) + sig*eps_t(it);
                    end
                end

                %% Estimate AR params
                xxx = X(1:K-1); yyy = X(2:K); 
                phi_hat = (xxx'*xxx)\xxx'*yyy;
                sig2_hat = sum( (yyy - xxx*phi_hat).^2)/(K-1);
                omega2_hat = sig2_hat/(1-phi_hat)^2;

                %% Estimate error process
                eps_hat = (X(K+t) - phi_hat*X(K+t-1))/sqrt(sig2_hat);

                %% Record AR params
                phi_save = phi_save + phi_hat;
                sig_save = sig_save + sqrt(sig2_hat);
                LRV_save = LRV_save + omega2_hat;

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
        
        twoSide_res = [twoSide_res;mean(detect_dummy)];
        oneSide_res = [oneSide_res;mean(detect_dummy2)];

    end
end

save(['Files/size_DGP',num2str(DGP),'_v01']);
%% Display
round(oneSide_res,2)