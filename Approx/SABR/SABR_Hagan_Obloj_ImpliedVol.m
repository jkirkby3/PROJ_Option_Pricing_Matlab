function IVs = SABR_Hagan_Obloj_ImpliedVol( Kvec, F0, nu, T, alpha, beta, rho )
% Calculates Implied volatility using asymptotic formulas from Obloj (correction to original SABR formulas)
% Kvec: strikes
% F0: initial forward price
% T: time to maturity of option
% (nu, alpha, beta, rho) are model parameters

IVs = zeros(length(Kvec),1);
for k = 1:length(Kvec)  % TODO: can easily vectorize this
    K = Kvec(k);
    
    z = (alpha / (nu*(1-beta)))*(F0^(1-beta) - K^(1-beta));
    chiz = log((sqrt(1 - 2*rho*z + z*z) + z - rho) / (1 - rho));

    FK = F0*K;  % NOTE: could reuse storage for K

    sig0 = alpha*log(F0/K)/chiz;
    sig1 = ((1-beta)*nu)^2 / 24 / FK^(1-beta) ...
           + 0.25*(rho*beta*alpha*nu) / FK^((1-beta)/2) ...
           + (2 - 3*rho^2) / 24 * alpha^2;

    IVs(k) = sig0 * (1.0 + sig1*T);

end

end

