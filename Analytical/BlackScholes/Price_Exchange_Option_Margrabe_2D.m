function price = Price_Exchange_Option_Margrabe_2D( S0_1, S0_2, T, rho, sigma_1, sigma_2, q_1, q_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: calcuates Price of 2D Exchange option, (S_1(T) - S_2(T)), using Magrabe Formula Under 2D Diffusion
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% S0_1    = initial asset price of first asset, e.g. S0_1 = 100
% S0_2    = initial asset price of second asset, e.g. S0_2 = 100
% T       = time to maturity in years (e.g. T=1)
% rho     = instantaneous correlation between S_1 and S_2
% sigma_1 = volatility of first asset (annualized), e.g. sigma = 0.2
% sigma_2 = volatility of second asset (annualized), e.g. sigma = 0.2
% q_1     = dividend yield of first asset, e.g. q_1 = 0.02
% q_2     = dividend yield of second asset, e.g. q_2 = 0.02
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig = sqrt(sigma_1^2 - 2*rho*sigma_1*sigma_2 + sigma_2^2);
d1 = (log(S0_1/S0_2) + 0.5*sig^2*T)/(sig*sqrt(T));
d2 = d1 - sig*sqrt(T);
price = exp(-q_1*T)*S0_1*normcdf(d1) - exp(-q_2*T)*S0_2*normcdf(d2);

end

