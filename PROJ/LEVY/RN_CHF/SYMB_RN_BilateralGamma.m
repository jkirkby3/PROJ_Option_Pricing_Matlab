function y = SYMB_RN_BilateralGamma( u, r, alpha_p, lam_p, alpha_m, lam_m)
% Risk-Neutral Levy Symbol for Bilateral Gamma model

zeta = -log((lam_p/(lam_p -1))^alpha_p*(lam_m/(lam_m +1))^alpha_m);
RNmu = r + zeta;

y = 1i*u*RNmu + ...
    log((lam_p./(lam_p -1i*u)).^alpha_p.*(lam_m./(lam_m + 1i*u)).^alpha_m);

end

