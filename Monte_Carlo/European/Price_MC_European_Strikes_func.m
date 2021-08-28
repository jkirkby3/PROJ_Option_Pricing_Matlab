function [prices, stdErrs] = Price_MC_European_Strikes_func(Spath, disc, call, Kvec )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Calculates European option prices for vector of strikes, given the simulatd paths 
% Returns: prices and standard errors for each of the supplied strikes
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% Spath = paths of underlying, dimension N_sim x M+1, where M = number of time steps (since includes S_0)
% disc = discount factor, e.g. exp(-r*T)
% call = 1 for call (else put)
% Kvec = strike vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prices = zeros(length(Kvec), 1);
stdErrs = zeros(length(Kvec), 1);
N_sim = size(Spath, 1);   % number of paths

for k = 1:length(Kvec)
    K = Kvec(k);
    if call ==1
        payoffs = max(0, Spath(:,end) - K);
    else
        payoffs = max(0, K - Spath(:,end));
    end
    prices(k) = disc*mean(payoffs);
    stdErrs(k) = disc*std(payoffs) / sqrt(N_sim);
end


end

