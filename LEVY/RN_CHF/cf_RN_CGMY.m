function y = cf_RN_CGMY(u,T,r,C,G,M,Y)
% CGMY RN CHF
% C,G,M>0 , Y <2
m = C*gamma(-Y)*((M-1)^Y - M^Y + (G+1)^Y - G^Y);  %Psi_s(-i)
y = C*T*gamma(-Y)*((M-1i*u).^Y - M^Y + (G+1i*u).^Y - G^Y); %T*Psi_s(u)
y = exp(1i*u*T*(r-m) + y);

end
