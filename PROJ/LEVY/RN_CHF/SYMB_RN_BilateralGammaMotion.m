function y = SYMB_RN_BilateralGammaMotion( u, r, alpha_p, lam_p, alpha_m, lam_m, sigma)
% Risk-Neutral Levy Symbol for Bilateral Gamma Motion model

sig2 = 0.5 * sigma^2;  % convexity correction for Brownian Motion component

zeta = -log((lam_p/(lam_p -1))^alpha_p*(lam_m/(lam_m +1))^alpha_m);
RNmu = r + zeta - sig2;

y = 1i*u*RNmu + ...
    log((lam_p./(lam_p -1i*u)).^alpha_p.*(lam_m./(lam_m + 1i*u)).^alpha_m);

% Add to the BG Model the Brownian Motion component
y = y - sig2 * u.^2;

end

