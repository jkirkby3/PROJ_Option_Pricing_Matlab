function [lx, r0, ux] = get_rate_grid_boundaries( model, modparam, t, gamma)
% Retrive lower and upper boundaries for variance grid, based on first two moments of variance process
% model : see below
% modparam : container with each of the model parameters required
% t: time to maturity used for the vairance grid
% gamma: grid mult param, rougly centered around +/- gamma standard deviations of the variance process

if model == 1 %CIR  (eta, theta, Rho, Sigmar, r0)
    eta = modparam.eta; theta = modparam.theta; Sigmar = modparam.Sigmar; r0 = modparam.r0; 

    mu_H = exp(-eta*t)*r0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmar^2/eta*r0*(exp(-eta*t)-exp(-2*eta*t)) +theta*Sigmar^2/(2*eta)*(1-2*exp(-eta*t)+exp(-2*eta*t));
    
elseif model == 2 %VASICEK ( in SV this is STEIN STEIN / same moments as SCOTT) ... (eta, theta, Rho, Sigmar, r0)
    eta = modparam.eta; theta = modparam.theta; Sigmar = modparam.Sigmar; r0 = modparam.r0; 
    
    mu_H = exp(-eta*t)*r0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmar^2/(2*eta)*(1 -exp(-2*eta*t));
    
elseif model == 3 % Merton
    lambda = modparam.lambda;  Sigmar = modparam.Sigmar; r0 = modparam.r0; 
    
    mu_H = r0 + lambda*t;
    sig2_H = Sigmar^2*t;
    
elseif model == 4 % Dothan
    a = modparam.a;  Sigmar = modparam.Sigmar; r0 = modparam.r0; 
    
    mu_H = r0*exp(a*t);
    sig2_H = r0^2*exp(2*a*t)*(exp(Sigmar^2*t) - 1);
    
end         

ux = mu_H + gamma*sqrt(sig2_H);  %variance grid upper bound
lx = mu_H - gamma*sqrt(sig2_H);

if model == 1 || model == 4
    lx = max(0.00001,lx);  %variance grid lower bound
end


end

