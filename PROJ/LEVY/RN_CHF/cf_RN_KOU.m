function y = cf_RN_KOU(u,T,r,sigma,lam,p_up,eta1,eta2)
% Kou Double Exponential RN CHF
% p_up = prob of upward jump
% lam = jump frequency

sig2 = .5*sigma^2;
temp1 = r - sig2 - lam*( (1-p_up)*eta2/(eta2+1) + p_up*eta1/(eta1-1) -1);
temp2 = -sig2*u.^2 + lam*( (1-p_up)*eta2./(eta2+1i*u) + p_up*eta1./(eta1-1i*u) - 1);
y = exp(T*1i*u*temp1 + T*temp2);

end
