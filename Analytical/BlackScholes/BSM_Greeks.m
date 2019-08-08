function [ Greeks ] = BSM_Greeks( G, S_0, sigma, r, q, T, K, call)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: calcuates prices and Greeks for Black Scholes Model (for strike or vector of strikes)
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% S_0   = initial asset price, e.g. S_0 = 100
% sigma = volatility (annualized), e.g. sigma = 0.2
% r     = interest rate, e.g. r = 0.05
% q     = dividend yeild, e.g. q = 0.02
% K     = strike (or vector of strikes)
% call  = 1 for call option (else put)
%
% G = Greek: 
%         0:Price
%         1:Delta, 2:Gamma, 3:Theta
%         4:Vega,  5:Rho,   6:Vanna
%         7:Vomma
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig = sigma;

eqT = exp(-q*T);
erT = exp(-r*T);
sqT = sqrt(T);

d1 = 1/(sig*sqT)*(log(S_0./K) + (r-q +sig^2/2)*T);
d2 = d1 - sig*sqT;

if G == 0 %Price
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
        Greeks = -eqT*S_0.*normpdf(d1)*sig/(2*sqT) - r*erT*K.*normcdf(d2)...
            + q*eqT*S_0.*normcdf(d1);
    else
        Greeks = -eqT*S_0.*normpdf(d1)*sig/(2*sqT) + r*erT*K.*normcdf(-d2)...
            -q*eqT*S_0.*normcdf(-d1);
    end
    
elseif G == 4 %Vegas
    Greeks = S_0.*normpdf(d1)*eqT*sqT;
    
elseif G == 5 %Rhos
    if call == 1
        Greeks = T*erT*K.*normcdf(d2);
    else
        Greeks = -T*erT*K.*normcdf(-d2);
    end
    
elseif G == 6 %Vannas
    Greeks = -eqT*normpdf(d1).*d2/sig;
    
elseif G == 7 %Vommas
    Greeks = eqT*sqT*S_0.*normpdf(d1).*d1.*d2/sig;
end

end

