function y = SYMB_RN_MJD(u,r,sigma, muj, sigmaj , lam)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
sig2 = .5*sigma^2; sigj2 = .5*sigmaj^2;
y = 1i*(r - sig2 - lam*(exp(muj + sigj2)-1))*u  -sig2*u.^2 + lam*(exp(1i*u*muj - sigj2*u.^2)-1);

end

