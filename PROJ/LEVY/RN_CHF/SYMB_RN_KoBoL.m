function y = SYMB_RN_KoBoL(u,r,c,lam_p,lam_m,nu)
% KoBoL RN Symbol - NOTE: params have been
% written in correspondence with CGMY, which is a subclass of KoBoL
C = c; M = lam_p; G = -lam_m; Y = nu;
m = C*gamma(-Y)*((M-1)^Y - M^Y + (G+1)^Y - G^Y);  %Psi_s(-i)
y = C*gamma(-Y)*((M-1i*u).^Y - M^Y + (G+1i*u).^Y - G^Y); %T*Psi_s(u)
y = 1i*u*(r-m) + y;

end
