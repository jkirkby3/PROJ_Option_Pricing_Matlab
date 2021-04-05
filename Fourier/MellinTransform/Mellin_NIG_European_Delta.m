function [ DELTA ] = Mellin_NIG_European_Delta( S_0, W, T, r, q, call, alpha, beta, delta, N1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N2 = N1;
N3 = N1;

f = W*exp(-r*T);

gam = sqrt(alpha^2 - beta^2);
k0 = log(S_0/W) + (r - q + delta*(sqrt(alpha^2 - (beta + 1)^2) - gam))*T;
adt = alpha*delta*T;
dta = 0.5*delta*T/alpha;

sum = 0;

if beta == 0
    % Symmetric Formula
    for n1 = 0:N1
        cons1 = k0^n1 / factorial(n1);
        for n2 = 1:N2
            term = besselk((n1-n2)/2+1, adt) / gamma((-n1+n2+1)/2)*(dta)^((-n1+n2)/2);
            sum = sum + cons1*term;
        end
    end
else
   % Asymmetric Formula
    
end

cons = f*alpha*exp(alpha*delta*T)/(S_0*sqrt(pi));
DELTA = cons*sum;


if call ~= 1  % delta of put using put-call parity
    DELTA = DELTA - exp(-q*T);
end

end


