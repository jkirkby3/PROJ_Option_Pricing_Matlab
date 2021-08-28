function [prices, stdErrs] = Price_MC_Asian_Strikes_func(Spath, call, disc, Kvec, M, mult)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Calculates arithmetic Asian option prices for vector of strikes, given the simulated paths 
% Returns: prices and standard errors for each of the supplied strikes
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% Spath = paths of underlying, dimension N_sim x M+1, where M = number of time steps (since includes S_0)
% call = 1 for call (else put)
% Kvec = strike vector
% M = number of monitoring points, e.g. 252 for "daily" monitoring
% mult = time partitioning multiplier to reduce bias (e.g. mult = 2 or 5)
% disc = discount factor, e.g. exp(-r*T), where T = Time (in years)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prices = zeros(length(Kvec),1);
stdErrs = zeros(length(Kvec),1);
N_sim = size(Spath,1);  % number of paths
M_mult = M*mult;  %time partitioning to reduce bias

Avg = zeros(N_sim,1); %Vector of payoffs, one per path
for n = 1:N_sim     
    Avg(n) = sum(Spath(n, 1:mult:M_mult+1));   %NOTE: mult is used to reduce bias
end
Avg = Avg /(M+1);

for k = 1:length(Kvec)
    K = Kvec(k);
    if call ==1
        payoffs = max(0, Avg - K);
    else
        payoffs = max(0, K - Avg);
    end
    
    prices(k) = disc*mean(payoffs);
    stdErrs(k) = disc*std(payoffs) / sqrt(N_sim);
end

end

