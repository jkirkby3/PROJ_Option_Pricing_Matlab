function y = cf_RN_BilateralGamma(u, r, T, alpha_p, lam_p, alpha_m, lam_m)
% Risk-Neutral Characterisitc Function for Bilateral Gamma model
zeta = -log((lam_p/(lam_p -1))^alpha_p*(lam_m/(lam_m +1))^alpha_m);
RNmu = r + zeta;

y = exp(T*(1i*u*RNmu))...
    .*(lam_p./(lam_p - 1i*u)).^(alpha_p*T)...
    .*(lam_m./(lam_m + 1i*u)).^(alpha_m*T);

end

