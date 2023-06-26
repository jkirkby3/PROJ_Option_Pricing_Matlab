function y = cf_RN_BilateralGammaMotion(u, r, T, alpha_p, lam_p, alpha_m, lam_m, sigma)
% Risk-Neutral Characterisitc Function for Bilateral Gamma model

sig2 = 0.5 * sigma^2;  % convexity correction for Brownian Motion component

zeta = -log((lam_p/(lam_p -1))^alpha_p*(lam_m/(lam_m +1))^alpha_m);
RNmu = r + zeta - sig2;

y = exp(T*(1i*u*RNmu - sig2 * u.^2))...
    .*(lam_p./(lam_p - 1i*u)).^(alpha_p*T)...
    .*(lam_m./(lam_m + 1i*u)).^(alpha_m*T);

end

