function [price, stdErr] = Price_MC_Var_Swaps_func(Spath, disc, M, mult)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Calculates Variance Swap prices (fair strike) for vector of strikes, given the simulated paths 
%           Uses Prices of Underlying for Control Variate
% Returns: prices and standard errors for each of the supplied strikes
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% Spath = paths of underlying, dimension N_sim x M+1, where M = number of time steps (since includes S_0)
% M = number of monitoring points, e.g. 252 for "daily" monitoring
% mult = time partitioning multiplier to reduce bias (e.g. mult = 2 or 5)
% disc = discount factor, e.g. exp(-r*T), where T = Time (in years)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N_sim = size(Spath,1);  % number of paths
M_mult = M*mult;  %time partitioning to reduce bias
S_0 = Spath(1,1);

RV = zeros(N_sim,1); %Vector of payoffs, one per path
for n = 1:N_sim  %NOTE: mult is used to reduce bias
    RV(n) = sum(log(Spath(n, 1+mult:mult:M_mult+1)./Spath(n, 1:mult:M_mult)).^2);
end

% Use S_T as Control Variate (CV) since E[S_T] = S_0/disc ... not great but provides a small reduction in std err

meanRV = mean(RV); 
price_NoCV  = disc*meanRV;

CV = Spath(:,M_mult + 1);

Ref = S_0/disc;
meanCV = mean(CV);
covXY = 1/(N_sim-1)*sum((CV - meanCV).*(RV - meanRV));
cstar = - covXY/var(CV);

price = price_NoCV + cstar*(meanCV - Ref);

stdErr = std(disc*RV + cstar*(CV - Ref))/sqrt(N_sim);

end

