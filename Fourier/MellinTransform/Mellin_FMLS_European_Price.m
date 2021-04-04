function [ price ] = Mellin_FMLS_European_Price( S_0, W, T, r, q, call, sigma, alpha, N1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N2 = N1;

F = W*exp(-r*T);
w_ = sigma^alpha / cos(pi*alpha/2);
k_ = log(S_0/W) + (r - q + w_)*T;
sum = 0;

wt = - w_ * T;

for n1 = 0:N1
    fn1 = factorial(n1);
    kn1 = k_^n1;
    for n2 = 1:N2
        c = (n1 - n2)/alpha;
        numer = kn1* wt^(-c);
        denom = fn1*gamma(1 - c);
        sum = sum + numer / denom;
    end
end

cons =  F/alpha;
price = cons*sum;

if call ~= 1  % price put using put-call parity
    price = price - (S_0*exp(-q*T) - W*exp(-r*T));
end

end

