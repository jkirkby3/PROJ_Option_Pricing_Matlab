 %Copyright (C) 2015 M.J. Ruijter

    %This file is part of BENCHOP.
    %BENCHOP is free software: you can redistribute it and/or modify
    %it under the terms of the GNU General Public License as published by
    %the Free Software Foundation, either version 3 of the License, or
    %(at your option) any later version.

    %BENCHOP is distributed in the hope that it will be useful,
    %but WITHOUT ANY WARRANTY; without even the implied warranty of
    %MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %GNU General Public License for more details.

    %You should have received a copy of the GNU General Public License
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.
function [U] = Price_Exchange_COS2D(S,T,r,sig1,sig2,rho,N,L)
% BENCHOP Problem 6: The Black-Scholes-Merton model for two underlying assets
% BSeuCallspread_COS computes the price for a European spread call option
%
% Input:    S       - Initial asset price   
%           K       - Strike price
%           T       - Terminal time  
%           r       - Risk-free interest rate
%           sig1    - Volatility process 1
%           sig2    - Volatility process 2
%           rho     - Correlation coefficient
%`          N       - number of cos expansion terms (e.g. 20)
%           L       - grid width parameter (e.g. 10)
%
% Output:   U       - Option value
%
% This MATLAB code has been written for the BENCHOP project and is based on 
% the COS methodes developed by F. Fang, C.W. Oosterlee, and M.J. Ruijter

if nargin < 7
    N = 19;
end
if nargin < 8
    L = 8;
end

% Parameters
sigma = [sig1 sig2];
rhoM = [1,rho;rho,1];
x = log(S);

% Parameters bivariate normal distribution
Sigma = (sigma'*sigma).*rhoM*T;
Mu = (r-0.5*sigma'.^2)*T;
% Interval [a,b]
a1 = log(100)+Mu(1)-L*sigma(1)*sqrt(T);
b1 = log(100)+Mu(1)+L*sigma(1)*sqrt(T);
a2 = log(100)+Mu(2)-L*sigma(2)*sqrt(T);
b2 = log(100)+Mu(2)+L*sigma(2)*sqrt(T);
a = min(a1,a2);
b = max(b1,b2);

% Number of Fourier cosine coefficients
N1 = N;
N2 = N1;
k1vec = [0:N1-1]';
k2vec = [0:N2-1];
omega1 = repmat(pi/(b-a)*[0:N2-1],N1,1);
omega(:,2) = omega1(1:end)';
omega(:,1) = repmat(pi/(b-a)*[0:N1-1]',N2,1);
omega_p = omega;
omega_m = omega;
omega_m(:,2) =-omega_m(:,2); 

% Fourier cosine coefficients payoff function
k1 = repmat(k1vec,1,N2);
k2 = repmat(k2vec,N1,1);
tempplus = (k1+k2)*pi;
tempmin = (k1-k2)*pi;
tempplusquad = tempplus.^2;
tempminquad = tempmin.^2;
bmina = b-a;
bminaquad = bmina^2;
piquad = pi^2;
k2quad = k2.^2;

Uk =-0.5*bmina^3./((tempplusquad+bminaquad).*(tempminquad+bminaquad).*k1.*k2*pi.*(piquad*k2quad+bminaquad)).*(...
        exp(b)*k1.*(tempplusquad+bminaquad).*(pi*k2.*tempmin+bminaquad).*sin(tempmin)...
        -exp(b)*k1.*(tempminquad+bminaquad).*(pi*k2.*tempplus+bminaquad).*sin(tempplus)...
        +exp(b)*pi*bmina*(2.*k2-k1).*k1.*(tempplusquad+bminaquad).*cos(tempmin)...     
        +exp(b)*pi*bmina*(2.*k2+k1).*k1.*(tempminquad+bminaquad).*cos(tempplus)...
    -2*exp(a)*k2.*(tempminquad.*tempplusquad.*sin(k1.*pi)+2*pi*bmina*k1.*(piquad*k2quad+bminaquad)) );
% for k1==k2
for i = 1:min(N1,N2);    Uk(i,i) = 0;   end
k2 = k1vec;
k2quad = k2.^2;
Ukdiag =-bminaquad./( (k2quad.*piquad+bminaquad).*k2*pi.*(4*piquad*k2quad+bminaquad) ).*...
    ( 0.5*exp(b)*bmina*(2*piquad*k2quad-bminaquad).*sin(2*pi*k2)...
   +1.5*exp(b)*pi*bminaquad.*k2.*cos(2*pi.*k2)...
-exp(a)*bmina*(4*piquad.*k2quad+bminaquad).*sin(pi*k2)...
-2*pi*k2.*( (-piquad*k2quad-0.25*bminaquad)*exp(b)+exp(a)*(piquad*k2quad+bminaquad) ) );
Uk = Uk+diag(Ukdiag,0);
% for (k1==0)&&(k2==0)
Uk(1,1) = (bmina+2)*exp(a)+exp(b)*(bmina-2);
% for k1==0
k2 = k2vec(2:end);
k2quad = k2.^2;
Uk(1,(2:end)) = bmina^3./( (piquad*k2quad+bminaquad).^2*pi.*k2 ).*...
    ( exp(b)*(pi*k2+bmina).*(-pi*k2+bmina).*sin(pi*k2)...
   +pi*k2.*(-2*exp(b)*bmina.*cos(pi*k2)+exp(a)*(piquad*k2quad+bmina*(bmina+2) ) ) );
% for k2==0
k1 = k1vec(2:end);
k1quad = k1.^2;
Uk(2:end,1) =  bmina./( (k1quad.*piquad+(bminaquad)).^2.*pi.*k1 ).*...
    ( (-(piquad*(-bmina+1)*k1quad+((-bmina+3)*bminaquad))*piquad*exp(b).*k1quad+exp(a)*(piquad*k1quad+(bminaquad)).^2 )...
.*sin(k1.*pi)+pi*bminaquad*k1.*((piquad*k1quad+((a-b+2)*(-b+a)))*exp(b).*cos(pi*k1)-2*exp(a)*- bmina) ) ;

Uk = 2/(b-a)*2/(b-a)*Uk;
Uk(:,1) = 0.5*Uk(:,1);    
Uk(1,:) = 0.5*Uk(1,:);  

% Characteristic function
cfplus = exp(1i*(Mu'*omega_p')-0.5*sum((omega_p').*(Sigma*omega_p'),1)); 
cfmin = exp(1i*(Mu'*omega_m')-0.5*sum((omega_m').*(Sigma*omega_m'),1));

% Fourier cosine coefficients density function
Recfplus = 1/2*real(repmat(cfplus,size(x,1),1).*exp(1i*((x-a)*omega_p')));
Recfmin = 1/2*real(repmat(cfmin,size(x,1),1).*exp(1i*((x-a)*omega_m')));
Recf = Recfplus+Recfmin;

% Option value 
U = exp(-r*T)*Uk(1:end)*Recf.';

end