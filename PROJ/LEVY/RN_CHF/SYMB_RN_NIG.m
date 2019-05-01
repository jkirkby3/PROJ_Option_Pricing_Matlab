function y = SYMB_RN_NIG(u,r,alph,bet,delt)
% Returns Risk Neutral SYMBOL
%   Detailed explanation goes here
asq = alph^2;
bsq = bet^2;
temp = sqrt(asq-bsq);
yy =  -delt*(sqrt(asq - (bet +1i*u).^2) - temp);  %Psi_s
RNmu = r + delt*(sqrt(asq - (bet+1)^2)-temp);
y = 1i*u*RNmu + yy;


end

