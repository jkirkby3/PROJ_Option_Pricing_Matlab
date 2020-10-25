function price = Price_Spread_Option_Kirk_2D(K, S0_1, S0_2, T, r, rho, sigma_1, sigma_2, q_1, q_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: calcuates Price of 2D Spread option, (S_1(T) - S_2(T) - K)+, using Kirks approximation
%        Note: when K = 0, this agrees with Magrabes exact formula
% Author: Justin Lars Kirkby
%
% Reference: Kirk, E. (1995): "Correlation in the energy markets," In Managing Energy Price
%               Risk (First Edition). London: Risk Publications and Enron, pp. 71-78.
%
% -----------------
% Params
% -----------------
% S0_1    = initial asset price of first asset, e.g. S0_1 = 100
% S0_2    = initial asset price of second asset, e.g. S0_2 = 100
% T       = time to maturity in years (e.g. T=1)
% r       = interest rate
% rho     = instantaneous correlation between S_1 and S_2
% sigma_1 = volatility of first asset (annualized), e.g. sigma = 0.2
% sigma_2 = volatility of second asset (annualized), e.g. sigma = 0.2
% q_1     = dividend yield of first asset, e.g. q_1 = 0.02
% q_2     = dividend yield of second asset, e.g. q_2 = 0.02
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F1 = S0_1*exp((r-q_1)*T);
F2 = S0_2*exp((r-q_2)*T);

f2k = F2 / (F2 + K);
sig = sqrt(sigma_1^2 - 2*rho*sigma_1*sigma_2*f2k + f2k^2*sigma_2^2);

d1 = (log(F1/(F1+K)) + 0.5*sig^2*T)/(sig*sqrt(T));
d2 = d1 - sig*sqrt(T);

price = exp(-r*T)*(F1*normcdf(d1) - (F2 + K)*normcdf(d2));

end

