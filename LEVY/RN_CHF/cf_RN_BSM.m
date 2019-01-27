function y = cf_RN_BSM( u, r, T, sigma )
% Return: risk neutral chf evaluated at u
% r = interest rate (or inputs r-q)
% sigma = volatility (per time unit?)
% T = time units til maturity

y = exp(T*(1i*(r-.5*sigma^2)*u - .5*sigma^2*u.^2));

end

