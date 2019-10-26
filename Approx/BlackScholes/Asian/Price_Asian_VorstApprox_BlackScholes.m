function price = Price_Asian_VorstApprox_BlackScholes(S_0, sigma, M, W, call, T, r, q, enforce_convention)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

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

