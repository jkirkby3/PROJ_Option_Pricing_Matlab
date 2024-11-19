function y = cf_RN_Heston_KDE(u,T,r,v_0,theta,kappa,sigma_v,rho, lam,p_up,eta1,eta2)
% HESTON + Kou Double Exponential Jumps Risk Neutral CHF  
% r = risk netural rate of interest (or pass r-q, interest minus div yield)
% T = time
% v_0 = spot variance
% theta = long term variance level
% kappa = rate of variance mean reversion
% sigma_v = volatility of variance
% rho = correlation between Brownian motions

% p_up = prob of upward jump
% lam = jump frequency

alpha = -.5*(u.^2 + u*1i);
beta = kappa - rho*sigma_v*u*1i;
omega2 = sigma_v^2;
gamma = .5*omega2;

D = sqrt(beta.^2 - 4.0*alpha.*gamma);

bD = beta - D;
eDt = exp(-D*T);

G = bD./(beta + D);
B = (bD./omega2).*((1.0-eDt)./(1.0 - G.*eDt));
psi = (1.0-G.*eDt)./(1.0-G);
A = ((kappa*theta)/(omega2))*(bD*T-2.0*log(psi));

temp1 = - lam*( (1-p_up)*eta2/(eta2+1) + p_up*eta1/(eta1-1) -1);
temp2 = lam*( (1-p_up)*eta2./(eta2+1i*u) + p_up*eta1./(eta1-1i*u) - 1);

%y = exp(A + B*v_0 + 1i*u*r*T);
y = exp(A + B*v_0 + T*1i*u*(r + temp1) + T*temp2);

end
