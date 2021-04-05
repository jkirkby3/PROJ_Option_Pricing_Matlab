function [ price ] = Mellin_VG_European_Price( S_0, W, T, r, q, call, sigma, theta, nu, N1, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European Options using Mellin Transform
% Models Supported: Variance Gamma Model, VG(sigma, theta, nu)
% Returns: price of contract
% Author: Justin Lars Kirkby/ Jean-Philippe Aguilar
%
% Reference: 1) "Some pricing tools for the Variance Gamma model", J-P Aguilar, 2020
%            2) "Pricing, risk and volatility in subordinated marketmodels", Aguilar, Kirkby, Korbel, 2020
%            3) "Closed-form option pricing in exponential Levy models", Aguilar and Kirkby, 2021
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
%
% sigma = param in VG model
% theta = param in VG model
% nu = param in VG model
%
% ----------------------
% Numerical Params 
% ----------------------
% N1  = maximum number summation terms in the series, will sum fewer terms
%       if error threshold (tol) is reached
% tol = desired error threshold of price (will stop adding terms once satisfied) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 11
    tol = 0;
end

if theta == 0
    price = Mellin_SymmetricVG_European_Price( S_0, W, T, r, q, call, sigma, nu, N1, tol);
    return;
end

F = W*exp(-r*T);

w_vg = log(1 - theta*nu - 0.5*sigma*sigma*nu) / nu;  % convexity correction
k_vg = log(S_0/W) + (r - q + w_vg)*T;

tau_v = T/nu - 0.5 + sqrt(2)*1e-13;   % add a little noise to prevent rational argument to gamma
sig_v = sigma*sqrt(nu/2);
theta_s = theta / (sigma*sigma);
%q_svt = (1/(sig_v*sig_v) + theta_s*theta_s)/4;
q_svt = 1/(2*sigma*sigma*nu) + (theta/(2*sigma*sigma))^2;

mult = F / (2^(2+2*tau_v) * sqrt(pi) * sig_v^(1+2*tau_v) * gamma(T/nu));
tol = tol / mult;

if k_vg < 0
    price = mult * price_minus( N1, tau_v, k_vg, q_svt, theta_s, tol);
    
elseif k_vg > 0
    p = mult * price_plus( N1, tau_v, k_vg, q_svt, theta_s, tol);
    price = S_0*exp(-q*T) - W*exp(-r*T) - p;
    
else
    price = 0; % TODO:
    % price = mult*price_zero( F, T, sigma, nu, N1);
end

if call ~= 1  % price put using put-call parity
    price = price - (S_0*exp(-q*T) - W*exp(-r*T));
end

end

function [a] = a_term(n1, n2, n3, tau_v, X, Y, Z)
% Computes a_{n1,n2,n3)(X,Y,Z)
    t2 = 2*tau_v;
    
    a = 2*gamma(-n1-tau_v)*gamma(-2*n1-n2-n3-t2-1)...
        *X^(2*n1+n2+n3+t2+1)...
        *Y^n1 /gamma(-2*n1-n2-t2);
    
    a = a + sqrt(pi)*gamma((-n1+n2+n3+t2+1)/2)...
        * pochhammer(-n1+n3+1,n2)...
        * X^n1 * Y^(-(-n1+n2+n3+t2+1)/2)...
        / ( 2^(-n1+n2+n3) * gamma(1+(-n1+n2+n3)/2) );
    
    a = a * Z^n2;
end

function [ sum ] = price_minus( N1, tau_v, k_vg, q_svt, theta_s, tol)
sum = 0;
last = 0;

for n1 = 0 : N1
    factn1 = factorial(n1);
    for n2 = 0 : N1
        m = (-1)^(n1 + n2) / (factn1 * factorial(n2));
        for n3 = 1 : N1
            p = a_term(n1, n2, n3, tau_v, -k_vg, q_svt, -theta_s);
            sum = sum + m*p;
        end
    end
    if n1 > 1 && abs(sum - last) < tol
        break;
    end
    last = sum;
end

% NOTE: the final mutiplier is done outside of this function
end

function [ sum ] = price_plus( N1, tau_v, k_vg, q_svt, theta_s, tol)
sum = 0;
last = 0;
for n1 = 0 : N1
    factn1 = factorial(n1);
    for n2 = 0 : N1
        cons = (factn1 * factorial(n2));
        for n3 = 1 : N1
            m = (-1)^(n1 + n2 + n3) / cons;
            p = a_term(n1, n2, n3, tau_v, k_vg, q_svt, theta_s);
            sum = sum + m*p;
        end
    end
    if n1 > 1 && abs(sum - last) < tol
        break;
    end
    last = sum;
end

% NOTE: the final mutiplier is done outside of this function
end

% function [p] = poch2(z, n)
%     p = gamma(z+n)/gamma(z);
% end

function p = pochhammer(a, n)
if (a == 0 && n <= 0) || (n == 0 && a > 0)
    p = 1;
elseif a == 0 && n > 0
    p = 0;
elseif a > 0
    if n == 1
        p = a;  % uses Gamma(a + 1) = a * Gamma(a)
    elseif n > 0
        p = prod(a:a + n - 1);
        % p = gamma(a + n)/gamma(a); 
    else
        p = inf; % TODO: what happens when a - n < 0
    end
else  
    p = neg_poch(a, n);
end
    
end

function p = neg_poch(m, n)
% Used for (-m)_n, m >= 1

m = -m;

if n > m
    p = 0;
else
    p = (-1)^n * factorial(m) / factorial(m - n);
end

end