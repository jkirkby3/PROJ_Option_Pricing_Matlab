function [ price ] = Mellin_SymmetricVG_European_Price( S_0, W, T, r, q, call, sigma, nu, N1, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European Options using Mellin Transform (Aguilar 2019)
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
% call  = 1 for call (else put)
% sigma = param in VG model
% nu = param in VG model
%
% ----------------------
% Numerical Params 
% ----------------------
% N1     = number summation terms in the series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 10
    tol = 0;
end

F = W*exp(-r*T);
k = log(S_0/F) - q*T;

theta = 0;
w_vg = log(1 - theta*nu - 0.5*sigma*sigma*nu) / nu;  % convexity correction
k_vg = k + w_vg*T;

if k_vg < 0
    price = price_minus( F, T, sigma, nu, N1, k_vg, tol);
    
elseif k_vg > 0
    p_minus = price_minus( F, T, -sigma, nu, N1, k_vg, tol);
    price = S_0*exp(-q*T) - W*exp(-r*T) - p_minus;
    
else
    price = price_zero( F, T, sigma, nu, N1, tol);
end

if call ~= 1  % price put using put-call parity
    price = price - (S_0*exp(-q*T) - W*exp(-r*T));
end

end


function [ sum ] = price_minus( F, T, sigma, nu, N1, k_vg, tol)
tau_vg = T/nu - 0.5 + sqrt(2)*1e-12;  % add a little noise to prevent rational argument to gamma
sigma_vg = sigma * sqrt(nu/2);

N2 = N1;

mult = F / (2*gamma(T/nu));
tol = tol/mult;

last = 0;
sum = 0;
for n1 = 0 : N1
    cons1 = (-1)^n1 / factorial(n1);
    cons2 = (-k_vg/sigma_vg)^n1;
    
    for n2 = 1 : N2
        t1 = gamma((-n1 + n2 + 1)/2 + tau_vg) / gamma((-n1 + n2)/2 + 1)...
            * cons2 * sigma_vg^n2;
        
        t2 = 2* gamma(-2*n1 - n2 - 1 - 2*tau_vg)/ gamma(-n1 + 0.5 - tau_vg) ...
            * (-k_vg/sigma_vg)^(2*n1 + 1 + 2*tau_vg) * (-k_vg)^n2;
        
        sum = sum + cons1*(t1 + t2);
    end
    if n1 > 1 && abs(sum - last) < tol
        break;
    end
    last = sum;
end

sum = sum * mult;
end


function [sum] = price_zero( F, T, sigma, nu, N1, tol)
tau_vg = T/nu - 0.5 + sqrt(2)*1e-12; % add a little noise to prevent rational argument to gamma
sigma_vg = sigma * sqrt(nu/2);

mult =  F / (2*gamma(T/nu));
tol = tol / mult;
last = 0;
sum = 0;
for n = 1 : N1
    t1 = sigma_vg^n * gamma((n+1)/2 + tau_vg) / gamma(n/2 + 1);
    sum = sum + t1;
    if n1 > 1 && abs(sum - last) < tol
        break;
    end
    last = sum;
end
sum = sum * mult;
end
