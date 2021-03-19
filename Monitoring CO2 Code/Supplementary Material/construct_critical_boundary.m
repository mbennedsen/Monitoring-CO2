%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "Designing a statistical procedure for monitoring global 
% carbon dioxide emissions" (2021) by Mikkel Bennedsen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code can be used to approximate the constant c_{T,alpha} for use in the boundary
% function, see the paper for details.
%
% INPUT:
% T: Monitoring horizon. We can approximate "indefinite" monitoring by
% setting T to a large value, e.g., T = 1000.
% alpha: Significance level (default: 5%, 10%, 32%).
% g_fct: Boundary function. This can be changed below to investigate
% alternative forms of the boundary function.
% R: Number of Monte Carlo replications.
%
% OUTPUT:
% c_tilde_oneSide: constant for use in one-sided hypothesis test.
% c_tilde_twoSide: constant for use in two-sided hypothesis test.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Mikkel Bennedsen, February 2021.
% Code can be freely used and distributed. Please cite Bennedsen (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
%% init
T = 30;% Number of years for monitoring

alpha = [0.05,0.10,0.32];   % Significance levels, alpha
R = 1e5;                    % Number of simulations

%% Boundary function (can be changed, see the Supplementary Matetial)
g_fct = @(t)( sqrt(t) ); %2
%% Calculate constant
c1 = nan(R,1);
c2 = nan(R,1);
for iR = 1:R
    W = cumsum(randn(T,1));
    c1(iR) = max( W./g_fct((1:T)') );
    
    c2(iR) = max( abs(W./g_fct((1:T)')) );
end
    
c_tilde_oneSide = [quantile(c1,1-alpha(1));
            quantile(c1,1-alpha(2));
            quantile(c1,1-alpha(3))];
        
c_tilde_twoSide = [quantile(c2,1-alpha(1));
            quantile(c2,1-alpha(2));
            quantile(c2,1-alpha(3))];

%% Save file
save(['Files/crit_val_T',num2str(T),'_v01'],'c_tilde_oneSide','c_tilde_twoSide');