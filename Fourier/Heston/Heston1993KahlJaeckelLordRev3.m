function [prices, alphas] = Heston1993KahlJaeckelLordRev3(PC, S,K,T,t,r,q,v0,theta,rho,kappa,sigma, alphas)
% Heston pricing function based on the implementation suggested by
% Roger Lord and Chrisitan Kahl in "Optimal Fourier inversion in 
% semi-analytical option pricing"
%
%
% Input: (PC till q can be vectorized)
%       PC: 1 for Calls, 2 for Puts
%       S: Spot
%       K: Strike
%       T: Maturity
%       t: start date
%       r: interest rate
%       q: dividend
%       v0: initial variance
%       theta: long run mean variance
%       kappa: mean reversion speed of  volatility
%       sigma: volatility of volatility
%       rho: correlation between returns volatility
%       alpha: alpha can be a vector supplied by the user, otherwise the
%       function attempts to find a payoff-dependent optimal alpha
%
%   Output: Price for each option, optionally generated alphas
%
%   Usage: Heston1993KahlJaeckelLordRev3(1, 100, 100, 20,0, 0.05, 0.0, 
%           0.00003, 0.00003,-0.3, 0.5, 0.0008)
%
%   Author: Jonathan Frei, 2015
% 

    % force column vector
    PC=PC(:);
    S=S(:);
    K=K(:);
    T=T(:);
    t=t(:);
    r=r(:);
    q=q(:);

    nos = numel(S);
    prices=NaN(nos,1);
    tau=T-t;
    mu=(r-q);
    F = S.*exp(mu.*tau);

    if(~exist('alphas','var'))
        alphas = NaN(numel(S),1);
    elseif(numel(alphas)==1)
        alphas = repmat(alphas,numel(S),1);
    end
    
    alpha0=0.75;
    
    for(ind=1:nos)
        if(isnan(alphas(ind)))
            try
                % using fzero here instead of fminsearch
               alphas(ind) = fzero( @(a) psi(a,K(ind), F(ind), kappa, theta, rho, sigma, tau(ind), v0), alpha0);
            catch
               alphas(ind) = alpha0;
            end
        end
        prices(ind) =  Ralpha(F(ind), K(ind), alphas(ind))+1/pi*integral(@(x) phi(x, K(ind), alphas(ind), F(ind), kappa, theta, rho, sigma, tau(ind), v0) , 0, Inf);
        if (PC(ind)==2)
            prices(ind) = prices(ind) + K(ind)*exp(-r(ind)*tau(ind))-S(ind)*exp(-q(ind)*tau(ind));
        end
    end

    
end


function p = psi(alpha, K, F, kappa, theta, rho, sigma, tau, v0)
    k = log(K);
    p = -alpha*k+0.5*log(phi(-(alpha+1)*1i, K, alpha, F, kappa, theta, rho, sigma, tau, v0)^2);
end

function r = Ralpha(F, K, alpha)
    r = F*(alpha<=0)-K*(alpha<=-1)-0.5*(F*(alpha==0)-K*(alpha==-1));
end

function y = phi(v, K, alpha, F, kappa, theta, rho, sigma, tau, v0)
    k = log(K);
    y = real(exp(-1i*(v-1i*alpha)*k).*( cf(v-1i*(alpha+1), F, kappa, theta, rho, sigma, tau, v0)./(-(v-1i*(alpha+1)).*(v-1i*alpha))));
end

function c = cf(u, F, kappa, theta, rho, sigma, tau, v0)
    f = log(F);
    c = exp(1i*u*f+ A(u, kappa, theta, rho, sigma, tau)+Bv(u, rho, sigma, kappa, tau)*v0);
end

function b = Bv(u, rho, sigma, kappa, tau)
    b = ((beta(u,rho,sigma,kappa)-D(u, rho, sigma, kappa)).*(1-exp(-D(u, rho, sigma, kappa)*tau)))./(sigma.^2*(1-G(u, rho, sigma, kappa).*exp(-D(u, rho, sigma, kappa)*tau)));
end

function a = A(u, kappa, theta, rho, sigma, tau)
    a = (kappa*theta*((beta(u,rho,sigma,kappa)-D(u, rho, sigma, kappa))*tau-2*log(phi2(u, rho, sigma, kappa, tau))))/sigma.^2;
end

function p = phi2(u, rho, sigma, kappa, tau)
    p = (G(u, rho, sigma, kappa).*exp(-D(u, rho, sigma, kappa)*tau)-1)./(G(u, rho, sigma, kappa)-1);
end

function g = G(u, rho, sigma, kappa)
    g = (beta(u,rho,sigma,kappa)-D(u, rho, sigma, kappa))./(beta(u,rho,sigma,kappa)+D(u, rho, sigma, kappa));
end

function d = D(u, rho, sigma, kappa)
    d = sqrt(beta(u,rho,sigma,kappa).^2-4*alphahat(u)*gamma(sigma));
end

function a = alphahat(u)
    a = -0.5*u.*(1i+u);
end

function b = beta(u,rho,sigma,kappa)
    b = kappa-rho*sigma*u*1i;
end

function y = gamma(sigma)
    y = 0.5*sigma.^2;
end