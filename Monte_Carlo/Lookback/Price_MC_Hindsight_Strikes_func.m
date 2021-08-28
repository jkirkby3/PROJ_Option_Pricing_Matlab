function [prices, stdErrs] = Price_MC_Hindsight_Strikes_func(Spath, call, strikes, M, mult, disc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Calculates Hindsight option prices for vector of strikes, given the simulatd paths 
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% Contract Types
% Fixed Strike (Hindsight) Put: (W - min{S_m: 0<=m<=M})^+
% Fixed Strike (Hindsight) Call: (max{S_m: 0<=m<=M} - W)^+
%
% Spath = paths of underlying, dimension N_sim x M+1, where M = number of time steps (since includes S_0)
% call = 1 for call (else put)
% Kvec = strike vector
% M = number of monitoring points, e.g. 252 for "daily" monitoring
% mult = time partitioning multiplier to reduce bias (e.g. mult = 2 or 5)
% S_0 = initial price
% disc = discount factor (e.g. exp(-r*T))
% T = Time (in years)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_sim = size(Spath,1);  % number of paths

M_mult = M*mult;  %time partitioning to reduce bias

prices = zeros(length(strikes),1);
stdErrs = zeros(length(strikes),1); %TODO: add code to return std errros

if call == 1
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

for k = 1:length(strikes)
    K = strikes(k);
    if call ==1  % Fixed Strike (Hindsight) Call: (max{S_m: 0<=m<=M} - W)^+
       payoffs  = max(0, curr_max - K);

    else  % Fixed Strike (Hindsight) Put: (W - min{S_m: 0<=m<=M})^+
        payoffs  = max(0, K - curr_min);
    end
    prices(k) = disc*mean(payoffs);
    stdErrs(k) = disc*std(payoffs) / sqrt(N_sim);
    
end

end

