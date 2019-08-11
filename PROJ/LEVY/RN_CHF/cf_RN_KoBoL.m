function y = cf_RN_KoBoL(u,T,r,c,lam_p,lam_m,nu)
% KoBoL RN CHF - NOTE: params have been
% written in correspondence with CGMY, which is a subclass of KoBoL
C = c; M = lam_p; G = -lam_m; Y = nu;
m = C*gamma(-Y)*((M-1)^Y - M^Y + (G+1)^Y - G^Y);  %Psi_s(-i)
y = C*T*gamma(-Y)*((M-1i*u).^Y - M^Y + (G+1i*u).^Y - G^Y); %T*Psi_s(u)
y = exp(1i*u*T*(r-m) + y);

end
