function [ price ] = Mellin_NegTempStable_European_Price( S_0, W, T, r, q, call, sigma, alpha, lambda, N1, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European Options using Mellin Transform
% Models Supported: Negative Tempered Stable
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
% lambda = param in model
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

N2 = N1;
N3 = N1;

F = W*exp(-r*T);
w_ = sigma^alpha / cos(alpha*pi/2);  % TODO: double check formula not same as appendix
wts_ = w_*((lambda + 1)^alpha - lambda^alpha);
c_ = - w_/gamma(-alpha);

k_ = log(S_0/W) + (r - q + wts_)*T;
sum = 0;
last = 0;
cons =  F*exp(lambda^alpha * w_*T) / alpha;
tol = tol / cons;

wt = - w_*T;
start_N3 = 1;  % Set to 1 for European, 0 for AON

for n1 = 0:N1
    fn1 = factorial(n1);
    for n2 = 0:N2
        fn2 = factorial(n2);
        for n3 = start_N3:N3
            g = pochhammer(1 - n1 + n3, n2);
            c = (n1 - n2 - n3)/alpha;
            term = g * lambda^n2 * k_^n1 * wt^(-c) / (fn1 * fn2 * gamma(1 - c) );
            sum = sum + term;
        end
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