function y = SYMB_RN_VG( u, r, sigma, nu, theta)
% Risk-Neutral Levy Symbol for Variance Gamma model
sig2 = .5*sigma^2; 
RNmu = r + log(1 - theta*nu - sig2*nu)/nu;
y = log(1 - 1i*theta*nu*u + sig2*nu*u.^2)/nu;
y = 1i*u*RNmu - y;

end

