function y = cf_RN_VG(u, r, T, sigma, nu, theta)
% Risk-Neutral Characterisitc Function for Variance Gamma model
sig2 = .5*sigma^2; 
RNmu = r + log(1 - theta*nu - sig2*nu)/nu;
y = log(1 - 1i*theta*nu*u + sig2*nu*u.^2)/nu;
y = exp(T*(1i*u*RNmu - y));

end

