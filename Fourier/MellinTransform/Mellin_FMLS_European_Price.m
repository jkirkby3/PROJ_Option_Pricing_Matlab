function [ price ] = Mellin_FMLS_European_Price( S_0, W, T, r, q, call, sigma, alpha, N1, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European Options using Mellin Transform
% Models Supported: Finite Moment Log Stable (FMLS)
% Returns: price of contract
% Author: Justin Lars Kirkby/ Jean-Philippe Aguilar
%
% Reference: 1) "Closed-form option pricing in exponential Levy models", Aguilar and Kirkby, 2021
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
% sigma = param in model
% alpha = param in model
%
% ----------------------
% Numerical Params 
% ----------------------
% N1  = maximum number summation terms in the series, will sum fewer terms
%       if error threshold (tol) is reached
% tol = desired error threshold of price (will stop adding terms once satisfied) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 10
    tol = 0;
end

N2 = N1;

F = W*exp(-r*T);
w_ = sigma^alpha / cos(pi*alpha/2);
k_ = log(S_0/W) + (r - q + w_)*T;
sum = 0; last = 0;

wt = - w_ * T;
cons =  F/alpha;
tol = tol / cons;

for n1 = 0:N1
    fn1 = factorial(n1);
    kn1 = k_^n1;
    for n2 = 1:N2
        c = (n1 - n2)/alpha;
        numer = kn1* wt^(-c);
        denom = fn1*gamma(1 - c);
        sum = sum + numer / denom;
    end
    if n1 > 1 && abs(sum - last) < tol
        break;
    end
    last = sum;
end


price = cons*sum;

if call ~= 1  % price put using put-call parity
    price = price - (S_0*exp(-q*T) - W*exp(-r*T));
end

end

