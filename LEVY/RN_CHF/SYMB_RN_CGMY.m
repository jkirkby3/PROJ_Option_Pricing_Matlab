function y = SYMB_RN_CGMY(u,r,C,G,M,Y)
% Returns Risk Neutral SYMBOL
% C,G,M>0 , Y <2
m = C*gamma(-Y)*((M-1)^Y - M^Y + (G+1)^Y - G^Y);  %Psi_s(-i)
y = C*gamma(-Y)*((M-1i*u).^Y - M^Y + (G+1i*u).^Y - G^Y); %T*Psi_s(u)
y = 1i*u*(r-m) + y;

end

