function [lx, v0, ux] = get_variance_grid_boundaries( model, modparam, t, gamma)
% Retrive lower and upper boundaries for variance grid, based on first two moments of variance process
% model : see below
% modparam : container with each of the model parameters required
% t: time to maturity used for the vairance grid
% gamma: grid mult param, rougly centered around +/- gamma standard deviations of the variance process

if model == 1 %HESTON  (eta, theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Sigmav = modparam.Sigmav; v0 = modparam.v0; 

    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmav^2/eta*v0*(exp(-eta*t)-exp(-2*eta*t)) +theta*Sigmav^2/(2*eta)*(1-2*exp(-eta*t)+exp(-2*eta*t));
    
elseif model == 2 || model == 6 %STEIN STEIN / SCOTT (eta, theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Sigmav = modparam.Sigmav; v0 = modparam.v0;
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmav^2/(2*eta)*(1 -exp(-2*eta*t));
    
elseif model == 3 % 3/2 MODEL
    %Transform to parameters that can be use in 4/2 model
    eta = modparam.eta*modparam.theta;
    theta = (modparam.eta + modparam.Sigmav^2)/eta;
    Sigmav = -modparam.Sigmav;
    v0 = 1/modparam.v0;
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmav^2/eta*v0*(exp(-eta*t)-exp(-2*eta*t)) +theta*Sigmav^2/(2*eta)*(1-exp(-eta*t)+exp(-2*eta*t));
    
elseif model == 4 % 4/2 MODEL (aa,bb, eta,theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Sigmav = modparam.Sigmav; v0 = modparam.v0;

    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmav^2/eta*v0*(exp(-eta*t)-exp(-2*eta*t)) +theta*Sigmav^2/(2*eta)*(1-exp(-eta*t)+exp(-2*eta*t));
    
elseif model == 5 % Hull White Model  (av, Rho, Sigmav, v0)
    Sigmav = modparam.Sigmav; v0 = modparam.v0; av = modparam.av;
    
    mu_H = v0*exp(av*t);
    sig2_H = v0^2*exp(2*av*t)*(exp(Sigmav^2*t)-1);
elseif model == 7 % alpha-Hypergeometric
%      %%% Commented Version is With the Second Formultion
%     Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; alphav = modparam.alphav; av = modparam.av; bv = modparam.bv;
%     Av     = 2*(av + Sigmav^2 - (v0)^(alphav/2)); 
%     mu_H   = v0*exp(Av*t);
%     sig2_H = v0^2*exp(2*av*t)*(exp((2*Sigmav)^2*t)-1);

    %%% Estimate Mean and Variance Using Stein Stein and First order approx
    Sigmav = modparam.Sigmav; v0 = modparam.v0; eta = modparam.eta; av = modparam.av; theta = modparam.theta;
    EtaBar = theta*av*exp(av*v0); ThetaBar = (eta - theta*exp(av*v0)*(1-av*v0))/EtaBar;
    mu_H = exp(-EtaBar*t)*v0 + ThetaBar*(1-exp(-EtaBar*t));
    %mu_H = v0;
    sig2_H = Sigmav^2/(2*EtaBar)*(1 -exp(-2*EtaBar*t));
end

if model == 8 %JACOBI
    lx = vmin ;
    ux = vmax;
else
    if model == 2 %For Stein-Stein, more sensitive
        lx = max(0.01,mu_H - gamma*sqrt(sig2_H));  %variance grid lower bound
    elseif model == 6 || model == 7   %For Scott/alpha-Hyper, allow lx negative
        lx = mu_H - gamma*sqrt(sig2_H);
    else
        lx = max(0.00001,mu_H - gamma*sqrt(sig2_H));  %variance grid lower bound
    end
    ux = mu_H + gamma*sqrt(sig2_H);  %variance grid upper bound
end


end

