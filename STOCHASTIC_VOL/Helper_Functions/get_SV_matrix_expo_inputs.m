function [ v1, v2, fv] = get_SV_matrix_expo_inputs( model,  modparam, psi_J, dt, v, dxi, r)
% 
if model == 1 %HESTON  (eta, theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav;

    c1 = (Rho*eta/Sigmav - .5);  c2 = (r - Rho*eta*theta/Sigmav);   c3 = .5*(1-Rho^2);
    v1 = dt*1i*(c1*v + c2 - psi_J(-1i));  %Note: we now have the compensated jump component
    v2 = dt*c3*v;
    fv = (1i*dxi*Rho/Sigmav)*v;
    
elseif model == 2  %STEIN STEIN (eta, theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav;
    
    c1 = Rho*eta/Sigmav - .5;  c2 = Rho*eta*theta/Sigmav;  c3 = r  - psi_J(-1i) - Rho*Sigmav/2; 
    vsq = v.^2;
    v1 = dt*1i*(c1*vsq - c2*v + c3);
    v2 = dt*.5*vsq*(1-Rho^2);
    fv = (1i*dxi*.5*Rho/Sigmav)*vsq;
    
elseif model == 3 || model == 4 % 3/2 MODEL or 4/2 model
    %Transform to parameters that can be use in 4/2 model
    if model == 3  % 3/2 (uses 4/2 model embedding)
        eta = modparam.eta*modparam.theta;
        theta = (modparam.eta + modparam.Sigmav^2)/eta;
        Sigmav = -modparam.Sigmav;
        % v0 = 1/modparam.v0;
        Rho = modparam.rho;
        aa = 0; bb = 1; %Now use the 4/2 Model embedding
    else
        eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav;
        aa = modparam.aa; bb = modparam.bb;
    end
        
    c1 = aa*Rho*eta/Sigmav - .5*aa^2; c2 = .5*(Rho*bb*Sigmav-bb^2) - bb*Rho*eta*theta/Sigmav;
    c3 = Rho*eta/Sigmav*(bb-aa*theta) + r - psi_J(-1i) - aa*bb; 
    sqrtv = sqrt(v);
    v1 = dt*1i*(c1*v + c2./v + c3);
    v2 = dt*.5*(1-Rho^2)*(aa*sqrtv + bb./sqrtv).^2;
    fv = (1i*dxi*Rho/Sigmav)*(aa*v + bb*log(v));
     
elseif model == 5 % Hull White Model  (av, Rho, Sigmav, v0)
    Rho = modparam.rho; Sigmav = modparam.Sigmav; av = modparam.av;
   
    c1 = .25*Rho*Sigmav - av*Rho/Sigmav;  c2 = .5;  c3 = r - psi_J(-1i);
    sqrtv = sqrt(v);
    v1 = dt*1i*(c1*sqrtv - c2*v + c3);
    v2 = dt*.5*(1-Rho^2)*v;
    fv = 1i*dxi*2*Rho/Sigmav*sqrtv;
    
elseif model == 6  % Scott
    eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav;

    c1v = Rho*(eta/Sigmav*(v - theta)-Sigmav/2);  %Note: this depends on v
    c2 = .5; c3 = r - psi_J(-1i);
    expv = exp(v); expv2 = expv.^2;
    v1 = dt*1i*(c1v.*expv - c2*expv2 + c3);
    v2 = dt*.5*(1-Rho^2)*expv2;
    fv = 1i*dxi*Rho/Sigmav*expv;
    
elseif model == 7 % alpha-Hypergeometric
%      %%% Commented Version is With the Second Formultion
%     Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; alphav = modparam.alphav; av = modparam.av; bv = modparam.bv;
%     
%     mu_func = @(v)2*(a+Sigmav^2)*v - 2*bv*v.^(1+alphav/2);
%     sig_func = @(v)2*Sigmav*v;
%
%     c1 = (Rho*Sigmav*.5 - Rho*(av + Sigmav)/Sigmav);
%     c2 = .5; c3 = Rho*bv/Sigmav; c4 = r - psi_J(-1i);
%     sqrtv = sqrt(v);
%     v1 = dt*1i*(c1*sqrtv  - c2*v + c3*v.^(1+alphav)/2 + c4 );
%     v2 = dt*.5*(1-Rho^2)*v;
%     fv = 1i*dxi*Rho/Sigmav*sqrtv;
  
    Rho = modparam.rho; Sigmav = modparam.Sigmav; eta = modparam.eta; av = modparam.av; theta = modparam.theta;

    c1 = Rho*theta/Sigmav;
    c2 = Rho*(eta/Sigmav + Sigmav/2); c3 = .5; c4 = r - psi_J(-1i);
    expv = exp(v); expv2 = expv.^2;
    v1 = dt*1i*(c1*exp((1+av)*v) - c2*expv - c3*expv2 + c4);
    v2 = dt*.5*(1-Rho^2)*expv2;
    fv = 1i*dxi*Rho/Sigmav*expv;
    
    
elseif model == 8 % JACOBI
    Rho = modparam.rho; Sigmav = modparam.Sigmav; eta = modparam.eta; theta = modparam.theta; vmin = modparam.vmin; vmax = modparam.vmax;
    %Qfunc = @(v) (v - vmin).*(vmax - v)/denomQ;
    %Qsqrt = @(v) sqrt((v - vmin).*(vmax - v)/denomQ);
    denomQ = (sqrt(vmax) - sqrt(vmin))^2;
    
    c1 = r - Rho*eta*theta/Sigmav - psi_J(-1i);
    c2 = Rho*eta/Sigmav - 0.5;
    v1 = dt*1i*(c1 + v*c2);
    v2 = dt*.5*(v - Rho^2/denomQ*(v - vmin).*(vmax - v));
    fv = (1i*dxi*Rho/Sigmav)*v;
    
%     c1 = (Rho*eta/Sigmav - .5);  c2 = (r - Rho*eta*theta/Sigmav);   c3 = .5*(1-Rho^2);
%     v1 = dt*1i*(c1*v + c2 - psi_J(-1i));  %Note: we now have the compensated jump component
%     v2 = dt*c3*v;
end

end

