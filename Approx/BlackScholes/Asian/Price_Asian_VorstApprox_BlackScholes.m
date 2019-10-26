function price = Price_Asian_VorstApprox_BlackScholes(S_0, sigma, M, W, call, T, r, q, enforce_convention)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Arithmetic Asian Options using Vorst Approximation Method
% Models Supported: Black-Scholes
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% ----------------------
% Contract/Model Params 
% ----------------------
% S_0   = initial stock price (e.g. 100)
% sigma = volatility of diffusion (e.g. 0.2)
% W     = strike  (e.g. 100)
% r     = interest rate (e.g. 0.05)
% q     = dividend yield (e.g. 0.05)
% T     = time remaining until maturity (in years, e.g. T=1)
% M     = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
% call  = 1 for call (else put)
% enforce_convention  =  true by default, enforces convention that S_0 is included in the average, else averaging starts at S_1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: we do an adjustment so that the formula is in terms of Avg(S_0,S_1,...,S_M) instead of Avg(S_1,...,S_M)

dt = T/M;

if nargin < 9
    enforce_convention = 1;
end

if enforce_convention
    W = (M+1)/M*W - S_0/M;  % Adjustment for avg convention
end

mu_G = log(S_0) + (r - q - 0.5*sigma^2)*(T + dt)/2;
sigma_G = sqrt(sigma^2*(dt + (T-dt)*(2*M - 1)/(6*M)));

if r - q == 0
    mult = M;
else
    mult = exp((r-q)*dt)*(1-exp((r-q)*M*dt))/(1-exp((r-q)*dt));
end

E_A = (S_0/M)*mult;
E_G = exp(mu_G + 0.5*sigma_G^2);
K = W - (E_A - E_G);  % Adjust strike based on difference between arithmetic and geometric

d1 = (mu_G - log(K) + sigma_G^2) / sigma_G;
d2 = d1 - sigma_G;

price = exp(-r*T)*(exp(mu_G + 0.5*sigma_G^2)*normcdf(d1) - K*normcdf(d2));

% Final adjustment so due to different averagin convention,  Avg(S_0,S_1,...,S_M) instead of Avg(S_1,...,S_M)
if enforce_convention
    price = price*(M/(M+1));
end

if call~=1  %Put Option
    if enforce_convention 
        if r - q == 0
            mult = M+1;
        else
            mult = (exp((r-q)*T*(1+1/M))-1)/(exp((r-q)*dt)-1);
        end
        price = price - S_0/(M+1)*exp(-r*T)*mult + W*exp(-r*T);  % NOTE: we use the ORIGINAL strike
    else
       price = price - S_0/(M)*exp(-r*T)*mult + W*exp(-r*T);  % NOTE: we use the ORIGINAL strike 
    end
    
end
end

