function [ mu_func,  sig_func] = get_SV_variance_grid_diffusion_funcs( model,  modparam)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if model == 1 %HESTON  (eta, theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Sigmav = modparam.Sigmav;
    mu_func  = @(v)eta*(theta-v);  %variance process drift
    sig_func = @(v)Sigmav*sqrt(v); %variance process vol
    
elseif model == 2  %STEIN STEIN (eta, theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Sigmav = modparam.Sigmav;
    mu_func  = @(v)eta*(theta-v);
    sig_func = @(v)Sigmav*(v>-100); %note the (v>-100) just ensures that sig_func returns a vector
    
    
elseif model == 3 || model == 4 % 3/2 MODEL or 4/2 model
    %Transform to parameters that can be use in 4/2 model
    if model == 3  % 3/2 (uses 4/2 model embedding)
        eta = modparam.eta*modparam.theta;
        theta = (modparam.eta + modparam.Sigmav^2)/eta;
        Sigmav = -modparam.Sigmav;
        % v0 = 1/modparam.v0;
    else
        eta = modparam.eta; theta = modparam.theta; Sigmav = modparam.Sigmav;
    end
        
    mu_func  = @(v)eta*(theta-v);
    sig_func = @(v)Sigmav*sqrt(v);
    
     
elseif model == 5 % Hull White Model  (av, Rho, Sigmav, v0)
    Sigmav = modparam.Sigmav; av = modparam.av;
    
    mu_func  = @(v)av*v;
    sig_func = @(v)Sigmav*v;

elseif model == 6  % Scott
    eta = modparam.eta; theta = modparam.theta; Sigmav = modparam.Sigmav;
    mu_func  = @(v)eta*(theta-v);
    sig_func = @(v)Sigmav*(v>-100); %note the (v>-100) just ensures that sig_func returns a vector

    
elseif model == 7 % alpha-Hypergeometric
%      %%% Commented Version is With the Second Formultion
%     Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; alphav = modparam.alphav; av = modparam.av; bv = modparam.bv;
%     
%     mu_func = @(v)2*(a+Sigmav^2)*v - 2*bv*v.^(1+alphav/2);
%     sig_func = @(v)2*Sigmav*v;

    Sigmav = modparam.Sigmav; eta = modparam.eta; av = modparam.av; theta = modparam.theta;
    mu_func = @(v)eta - theta*exp(av*v);
    sig_func = @(v)Sigmav*(v>-100);
    
    
elseif model == 8 % JACOBI
    Sigmav = modparam.Sigmav; eta = modparam.eta; theta = modparam.theta; vmin = modparam.vmin; vmax = modparam.vmax;
    %Qfunc = @(v) (v - vmin).*(vmax - v)/denomQ;
    %Qsqrt = @(v) sqrt((v - vmin).*(vmax - v)/denomQ);
    denomQ = (sqrt(vmax) - sqrt(vmin))^2;
    mu_func  = @(u) eta*(theta - u);
    sig_func = @(u) Sigmav/sqrt(denomQ)*sqrt((u - vmin).*(vmax - u));
    
end

end

