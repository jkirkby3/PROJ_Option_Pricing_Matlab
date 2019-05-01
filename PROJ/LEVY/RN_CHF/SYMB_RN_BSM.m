function y = SYMB_RN_BSM( u, r, sigma )
% Return: risk neutral SYMBOL evaluated at u
% r = interest rate (or inputs r-q)
% sigma = volatility (per time unit?)

y = 1i*(r-.5*sigma^2)*u - .5*sigma^2*u.^2;

end

