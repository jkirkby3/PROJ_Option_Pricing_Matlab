function [price, stdErr] = Price_MC_Lookback_func(Spath, call, M, mult, disc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Calculates Looback option price, given the simulatd paths 
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% Contract Types
% Floating Strike (Lookback) Put:  max{S_m: 0<=m<=M} - S_T
% Floating Strike (Lookback) Call: S_T - min{S_m: 0<=m<=M}
%
% Spath = paths of underlying, dimension N_sim x M+1, where M = number of time steps (since includes S_0)
% call = 1 for call (else put)
% M = number of monitoring points, e.g. 252 for "daily" monitoring
% mult = time partitioning multiplier to reduce bias (e.g. mult = 2 or 5)
% S_0 = initial price
% disc = discount factor (e.g. exp(-r*T))
% T = Time (in years)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_sim = size(Spath,1);  % number of paths

M_mult = M*mult;  %time partitioning to reduce bias


if call ~= 1
    curr_max = zeros(N_sim, 1);
    for n = 1:N_sim
        curr_max(n) = max(Spath(n, 1:mult:M_mult+1));
    end
else  % find_min
    curr_min = zeros(N_sim, 1);
    for n = 1:N_sim
        curr_min(n) = min(Spath(n, 1:mult:M_mult+1));
    end
end


if call ==1  % Floating Strike (Lookback) Call: S_T - min{S_m: 0<=m<=M}
   payoffs  = Spath(:, M_mult+1) - curr_min;

else  % Floating Strike (Lookback) Put:  max{S_m: 0<=m<=M} - S_T
    payoffs  = curr_max - Spath(:,M_mult+1);
end
price = disc*mean(payoffs);
stdErr = disc*std(payoffs) / sqrt(N_sim);

end

