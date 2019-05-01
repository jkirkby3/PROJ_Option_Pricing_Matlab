function [KdHW, KcHW] = hullwhitefairstrike(r, V0, sigma, mu, T, rho, n)
% From C. Bernard and Z. Cui, Prices and Asymptotics for Discrete Variance
% Swaps, Applied Mathematical Finance 2014, 21(2), 140-173

%%%%%%%%%%%%%%% for Hull-White - Proposition 4.1 %%%%%

KcHW=continuousHW(mu,T,V0);
KdHW=discreteHW(mu,T,rho,sigma,V0,r,n);                                 
  
function res=continuousHW(mu,T,V0)
    res= V0 * (exp(T * mu) - 1) / T / mu;
end

function res=discreteHW(mu,T,rho,sigma,V0,r,n)
    Kc=continuousHW(mu,T,V0);
    res=( (r^2*T)/n+ ( 1-r*T/n)*Kc -V0^2*(exp((2*mu+sigma^2)*T)-1)*(exp(mu*T/n)-1)/(2*T*mu*(mu+sigma^2)*(exp((2*mu+sigma^2)*T/n)-1))+V0^2*(exp((2*mu+sigma^2)*T)-1)/(2*T*(2*mu+sigma^2)*(mu+sigma^2))...
        +8*rho*( exp(3*(4*mu+sigma^2)*T/8)-1 )*V0^(3/2)*sigma*(exp(mu*T/n)-1)/(mu*T*(4*mu+3*sigma^2)*(exp(3*(4*mu+sigma^2)*T/(8*n))-1))-64*rho*(exp(3*(4*mu+sigma^2)*T/8)-1)*V0^(3/2)*sigma/(3*T*(4*mu+sigma^2)*(4*mu+3*sigma^2)));
end


end