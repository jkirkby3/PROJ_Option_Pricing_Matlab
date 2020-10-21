function [ price ] = Mellin_SymmetricVG_European_Price( S_0, W, T, r, q, call, sigma, nu, N1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European Options using Mellin Transform (Aguilar 2019)
% Models Supported: Symmetric Variance Gamma Model, VG(sigma, 0, nu)
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% Reference: "Some pricing tools for the Variance Gamma model", J-P Aguilar, 2020
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


F = W*exp(-r*T);
k = log(S_0/F) - q*T;

theta = 0;
w_vg = log(1 - theta*nu - 0.5*sigma*sigma*nu) / nu;  % convexity correction
k_vg = k + w_vg*T;

if k_vg < 0
    price = price_minus( F, T, sigma, nu, N1, k_vg);
    
elseif k_vg > 0
    p_minus = price_minus( F, T, -sigma, nu, N1, k_vg);
    price = S_0*exp(-q*T) - W*exp(-r*T) - p_minus;
    
else
    price = price_zero( F, T, sigma, nu, N1);
end

if call ~= 1  % price put using put-call parity
    price = price - (S_0*exp(-q*T) - W*exp(-r*T));
end

end


function [ price ] = price_minus( F, T, sigma, nu, N1, k_vg)
tau_vg = T/nu - 0.5;
sigma_vg = sigma * sqrt(nu/2);

N2 = N1;

price = 0;
for n1 = 0 : N1
    cons1 = (-1)^n1 / factorial(n1);
    cons2 = (-k_vg/sigma_vg)^n1;
    
    for n2 = 1 : N2
        t1 = gamma((-n1 + n2 + 1)/2 + tau_vg) / gamma((-n1 + n2)/2 + 1)...
            * cons2 * sigma_vg^n2;
        
        t2 = 2* gamma(-2*n1 - n2 - 1 - 2*tau_vg)/ gamma(-n1 + 0.5 - tau_vg) ...
            * (-k_vg/sigma_vg)^(2*n1 + 1 + 2*tau_vg) * (-k_vg)^n2;
        
        price = price + cons1*(t1 + t2);
    end
end

price = price * F / (2*gamma(T/nu));
end


function [price] = price_zero( F, T, sigma, nu, N1)
tau_vg = T/nu - 0.5;
sigma_vg = sigma * sqrt(nu/2);

price = 0;
for n = 1 : N1
    t1 = sigma_vg^n * gamma((n+1)/2 + tau_vg) / gamma(n/2 + 1);
    price = price + t1;
end
price = price * F / (2*gamma(T/nu));
end
