function y = SYMB_RN_Kou( u, r, sigma,lam,p_up,eta1,eta2)
% Kou Double Exponential RN symbol
% p_up = prob of upward jump

sig2 = .5*sigma^2;
temp1 = r - sig2 - lam*( (1-p_up)*eta2/(eta2+1) + p_up*eta1/(eta1-1) -1);
temp2 = -sig2*u.^2 + lam*( (1-p_up)*eta2./(eta2+1i*u) + p_up*eta1./(eta1-1i*u) - 1);
y = 1i*u*temp1 + temp2;

end

