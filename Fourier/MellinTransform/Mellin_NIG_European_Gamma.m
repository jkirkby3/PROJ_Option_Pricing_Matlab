function [ GAMMA ] = Mellin_NIG_European_Gamma( S_0, W, T, r, q, alpha, beta, delta, N1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

f = W*exp(-r*T);

gam = sqrt(alpha^2 - beta^2);
k0 = log(S_0/W) + (r - q + delta*(sqrt(alpha^2 - (beta + 1)^2) - gam))*T;
adt = alpha*delta*T;
dta = 0.5*delta*T/alpha;

sum = 0;

if beta == 0
    % Symmetric Formula
    for n = 0:N1
        cons1 = k0^n / factorial(n);
        term = besselk(n/2 + 1, adt) / gamma((-n+1)/2) * (dta)^(-n/2);
        sum = sum + cons1*term;
    end
else
   % Asymmetric Formula
    
end

cons = f*alpha*exp(alpha*delta*T)/(S_0*S_0*sqrt(pi));
GAMMA = cons*sum;


end


