function [ Greeks ] = BSM_Greeks( G, S_0, sigma, r, q, T, K, call)
% G = Greek: 
%         0:Price
%         1:Delta, 2:Gamma, 3:Theta
%         4:Vega,  5:Rho,   6:Vanna
%         7:Vomma
%  
sig = sigma;

eqT = exp(-q*T);
erT = exp(-r*T);
sqT = sqrt(T);

d1 = 1/(sig*sqT)*(log(S_0./K) + (r-q +sig^2/2)*T);
d2 = d1 - sig*sqT;


Npdf = @(x)exp(-x.^2/2)/sqrt(2*pi);

if G == 0 %PRICE
    if call == 1
        Greeks = eqT*S_0.*normcdf(d1) - erT*K.*normcdf(d2);
    else
        Greeks = erT*K.*normcdf(-d2) - eqT*S_0.*normcdf(-d1);
    end
    
elseif G==1 %Deltas
    if call == 1
        Greeks = eqT*normcdf(d1);
    else
        Greeks = - eqT*normcdf(-d1);
    end
elseif G==2 %Gammas
    Greeks = exp(-.5*d1.^2)./(sqrt(2*pi*T)*sig*S_0);
elseif G==3 %Thetas
    if call == 1
        Greeks = -eqT*S_0.*Npdf(d1)*sig/(2*sqT) - r*erT*K.*normcdf(d2)...
            + q*eqT*S_0.*normcdf(d1);
    else
        Greeks = -eqT*S_0.*Npdf(d1)*sig/(2*sqT) + r*erT*K.*normcdf(-d2)...
            -q*eqT*S_0.*normcdf(-d1);
    end
elseif G == 4 %Vegas
    Greeks = S_0.*Npdf(d1)*eqT*sqT;
elseif G == 5 %Rhos
    if call == 1
        Greeks = T*erT*K.*normcdf(d2);
    else
        Greeks = -T*erT*K.*normcdf(-d2);
    end
elseif G == 6 %Vannas
    Greeks = -eqT*Npdf(d1).*d2/sig;
elseif G == 7 %Vommas
    Greeks = eqT*sqT*S_0.*Npdf(d1).*d1.*d2/sig;
end
end

