function y = cf_RN_MJD( u, r, T, sigma, muj, sigmaj , lam)
% Return: risk neutral chf evaluated at u
% r = interest rate
% sigma = volatility (per time unit?)
% T = time units til maturity
% Jump ~ N(muj,sigmaj^2), i.e jumps in log return are normal
% lam = E[#jumps per unit of T]

sig2 = .5*sigma^2;
sigj2 = .5*sigmaj^2;
y = exp( 1i*(r - sig2 - lam*(exp(muj + sigj2)-1))*u*T  ...
        -sig2*u.^2*T + lam*T*(exp(1i*u*muj - sigj2*u.^2)-1));

end

