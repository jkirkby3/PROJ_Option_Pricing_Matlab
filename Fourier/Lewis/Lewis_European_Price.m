function [ price ] = Lewis_European_Price(S_0, W, rnCHF, T, r, q, call, max_u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European Options using the method of Lewis (2001)
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
% Author: Justin Lars Kirkby
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
% rnCHF = risk netural characteristic function (function handle with single argument)
% 
% ----------------------
% Numerical Params 
% ----------------------
% max_u = upper limit of fourier integral for numerical integration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 8
    max_u = 400;
end

k = log(S_0./W);

% Integrate the complex integral
grand = @(u)  real(exp(1i*u*k) .* rnCHF(u - 0.5*1i) ./ (u.*u + 0.25));
z = integral(grand, 0, max_u);

% Call Price
price = S_0*exp(-q*T) - (1/pi)*sqrt(S_0*W) * exp(-r*T) * z ;

if call ~= 1  % price put using put-call parity
    price = price - (S_0*exp(-q*T) - W*exp(-r*T));
end

end
