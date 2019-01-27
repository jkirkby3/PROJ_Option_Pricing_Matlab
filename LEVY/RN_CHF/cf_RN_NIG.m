function y = cf_RN_NIG( u,r,T,alpha,beta,delta)
%
%  
asq = alpha^2;
bsq = beta^2;
temp = sqrt(asq-bsq);
y = -delta*(sqrt(asq - (beta +1i*u).^2) - temp);  %Psi_s
RNmu = r + delta*(sqrt(asq - (beta+1)^2)-temp);
y = exp(T*(1i*u*RNmu + y));

end

