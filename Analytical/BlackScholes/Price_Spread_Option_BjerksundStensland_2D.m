function price = Price_Spread_Option_BjerksundStensland_2D(K, S0_1, S0_2, T, r, rho, sigma_1, sigma_2, q_1, q_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: calcuates Price of 2D Spread option, (S_1(T) - S_2(T) - K)+, using Bjerksund & Stensland approximation
%        Note: when K = 0, this agrees with Magrabes exact formula
% Author: Justin Lars Kirkby
%
% Reference: Bjerksund, P. and Stensland, G. (2006): "Closed form spread option valuation"
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

a = F2 + K;
b = F2 / a;
rhosigs = rho*sigma_1*sigma_2;
sig = sqrt(sigma_1^2 - 2*rhosigs*b + b^2*sigma_2^2);
sigst = sig*sqrt(T);


d1 = (log(F1/a) + (0.5*sigma_1^2 - b*rhosigs + 0.5*b^2*sigma_2^2)*T) / sigst;
d2 =  (log(F1/a) + (-0.5*sigma_1^2 + rhosigs + (0.5*b^2 - b)*sigma_2^2)*T) / sigst;
d3 = (log(F1/a) + (-0.5*sigma_1^2 + 0.5*b^2*sigma_2^2)*T) / sigst;

price = exp(-r*T)*(F1*normcdf(d1) - F2*normcdf(d2) - K*normcdf(d3));

end

