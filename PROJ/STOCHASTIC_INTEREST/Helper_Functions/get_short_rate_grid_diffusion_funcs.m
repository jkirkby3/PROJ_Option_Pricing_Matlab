function [ mu_func,  sig_func, zeta_func, varrho_func, nu, beta, is_affine_varrho] = get_short_rate_grid_diffusion_funcs( model,  modparam)
% Get the drift and diffusion function for one-factor short rate model
%  	dr_t = mu_func(r)*dt + sig_func(r)*dW_t
%
rho = modparam.rho;
Sigma_s = modparam.Sigma_s;
zeta_func = -1;
varrho_func = -1;

% Affine function params
nu = 0;
beta = 0;
is_affine_varrho = 0;

if model == 1 % CIR  (eta, theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Sigmar = modparam.Sigmar;
    mu_func  = @(v)eta*(theta-v);  %variance process drift
    sig_func = @(v)Sigmar*sqrt(v); %variance process vol
    
    c = rho*Sigma_s/Sigmar;
    
    zeta_func = @(r0, r1) 2*c * (sqrt(r1) - sqrt(r0));
    varrho_func = @(r) c*(eta*(theta - r) - Sigmar^2/4)./sqrt(r);
    
elseif model == 2  % Vasicek / Hull-White  (Note: in SV this is Stein-Stein)
    is_affine_varrho = 1;
    eta = modparam.eta; theta = modparam.theta; Sigmar = modparam.Sigmar; 
    
    c = rho*Sigma_s/Sigmar;
    
    mu_func  = @(v)eta*(theta-v);
    sig_func = @(v)Sigmar*(v>-1000); %note the (v>-100) just ensures that sig_func returns a vector
    zeta_func = @(r0, r1) c * (r1 - r0);
    varrho_func = @(r) c*eta*(theta - r);
    
    nu = c * eta * theta;
    beta = -c * eta;

elseif model == 3  % Merton
    is_affine_varrho = 1;
    lambda = modparam.lambda; Sigmar = modparam.Sigmar; 
    
    c = rho*Sigma_s/Sigmar;
    nu = lambda*c;
    beta = 0;
    mu_func  = @(v)lambda*(v>-1000);
    sig_func = @(v)Sigmar*(v>-1000);
    zeta_func = @(r0, r1) c * (r1 - r0);
    varrho_func = @(r) nu*(r>-1000);
    
elseif model == 4  % Dothan
    is_affine_varrho = 1;
    a = modparam.a; Sigmar = modparam.Sigmar; 
    
    c = rho*Sigma_s/Sigmar;
    nu = c*(a - Sigma_s^2/2);
    beta = 0;
    mu_func  = @(v)a*v;
    sig_func = @(v)Sigmar*v;
    zeta_func = @(r0, r1) c * log(r1/r0);
    varrho_func = @(r) nu*(r>-1000);
    
end

