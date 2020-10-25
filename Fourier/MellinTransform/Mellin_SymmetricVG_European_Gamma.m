function [ gam ] = Mellin_SymmetricVG_European_Gamma( S_0, W, T, r, q, sigma, nu, N1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Gamma Function for European Options using Mellin Transform (Aguilar, Kirkby, Korbel 2020)
% Models Supported: Symmetric Variance Gamma Model, VG(sigma, 0, nu)
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% Reference: 1) "Some pricing tools for the Variance Gamma model", J-P Aguilar, 2020
%            2) "Pricing, risk and volatility in subordinated marketmodels", Aguilar, Kirkby, Korbel, 2020
%
% ----------------------
% Contract/Model Params 
% ----------------------
% S_0 = initial stock price (e.g. 100)
% W   = strike  (e.g. 100)
% r   = interest rate (e.g. 0.05)
% q   = dividend yield (e.g. 0.05)
% T   = time remaining until maturity (in years, e.g. T=1)
% sigma = param in VG model
% nu = param in VG model
%
% ----------------------
% Numerical Params 
% ----------------------
% N1     = number summation terms in the series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = T + 1e-7*sqrt(2);  % add a little noise to protect against rational singularities

f = W*exp(-r*T);
k = log(S_0/f) - q*T;

theta = 0;
w_vg = log(1 - theta*nu - 0.5*sigma*sigma*nu) / nu;  % convexity correction
k_vg = k + w_vg*T;

if k_vg < 0
    gam = gam_minus( S_0, f, T, sigma, nu, N1, k_vg);
    
elseif k_vg > 0
    gam = -gam_minus( S_0, f, T, -sigma, nu, N1, k_vg);
    
else
    gam = gam_zero( S_0, f, T, sigma, nu);
end

end


function [ gam ] = gam_minus( S, f, T, sigma, nu, N1, k_vg)
tau_vg = T/nu - 0.5 + sqrt(2)*1e-12;  % add a little noise to prevent rational argument to gamma
sigma_vg = sigma * sqrt(nu/2);

gam = 0;
for n = 0 : N1

    t1 = gamma(-n/2 + tau_vg) / gamma((-n + 1)/2) * (-k_vg/sigma_vg)^(n) ;
    t2 = 2*gamma(-2*n - 2*tau_vg)/gamma(-n + 0.5 - tau_vg) * (-k_vg/sigma_vg)^(2*n + 2*tau_vg);
    
    cons1 = (-1)^n / factorial(n);
    gam = gam + cons1*(t1 + t2);
end

gam = gam * f / (2 * S*S * sigma_vg * gamma(T/nu));
end


function [gam] = gam_zero( S, f, T, sigma, nu)
tau_vg = T/nu - 0.5 + sqrt(2)*1e-12; % add a little noise to prevent rational argument to gamma
sigma_vg = sigma * sqrt(nu/2);

gam = f / (2*sqrt(pi) * S*S * sigma_vg * gamma(T/nu)) * gamma(tau_vg - 0.5)/gamma(T/nu);
end
