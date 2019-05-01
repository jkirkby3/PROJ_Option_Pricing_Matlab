function alph = GetAlph_DisreteVariance( c2Jump, c4Jump, model, modparam, T, L1 )
% Gets alph parameter (truncated density width param) to apply PROJ method to ChF of realized variance
% c2Jump - second cumulant of jump component
% c4Jump - fourth cumulant of jump component
% model - SV model id number (see below)
% modparm - container with each param needed for given model
% T - contract maturity
% L1 - truncation gridwidth multiplication param

minalph = .5;  % minimun allowable alpha

t = T/2;

if model == 1 %HESTON  (eta, theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; 
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    
    mfunc = @(v)sqrt(v);
    
elseif model == 2 %STEIN STEIN 
    eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0;
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    
     mfunc = @(v) v;
elseif model == 6 %%SCOTT
    eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0;
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    
     mfunc = @(v)exp(v);
     
elseif model == 3 % 3/2 MODEL
    %Transform to parameters that can be use in 4/2 model
    eta = modparam.eta*modparam.theta;
    theta = (modparam.eta + modparam.Sigmav^2)/eta;
    Sigmav = -modparam.Sigmav;
    v0 = 1/modparam.v0;
    Rho = modparam.rho;
    aa = 0; bb = 1; %Now use the 4/2 Model:
     
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    mfunc = @(v)1./sqrt(v);
    
elseif model == 4 % 4/2 MODEL (aa,bb, eta,theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0;
    aa = modparam.aa; bb = modparam.bb;
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    mfunc = @(v) aa*sqrt(v) + bb/sqrt(v);
    
elseif model == 5 % Hull White Model  (av, Rho, Sigmav, v0)
    Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; av = modparam.av;
    
    mu_H = v0*exp(av*t);
    mfunc = @(v)sqrt(v);
    
elseif model == 7 % alpha-Hypergeometric
    Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; eta = modparam.eta; av = modparam.av; theta = modparam.theta;
    
    %%% Estimate Mean and Variance Using Stein Stein and First order approx
    EtaBar = theta*av*exp(av*v0); ThetaBar = (eta - theta*exp(av*v0)*(1-av*v0))/EtaBar;
    mu_H = exp(-EtaBar*t)*v0 + ThetaBar*(1-exp(-EtaBar*t));
    mfunc = @(v)exp(v);
end

c2 = c2Jump + mfunc(mu_H)^2;
c4 = c4Jump;
alph = L1*sqrt(c2*t + sqrt(c4*T));
alph = max(alph, minalph);
end

