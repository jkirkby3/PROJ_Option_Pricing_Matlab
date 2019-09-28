function Spath = Simulate_Heston_Euler_Schemes( N_sim, M, T, S_0, r, q, SVModelParams, scheme)
% Simulates Paths of Heston Model under various Euler Schemes
% See Lord et al (2010) for details
% N_sim = # paths
% M = #time steps on [0,T], ie dt =T/M   
% Note: returns paths of dimension (N_sim,M+1), since they include S_0
%
% scheme:     
%        1 = Absorbption
%        2 = Reflection
%        3 = Higham and Mao
%        4 = Partial Truncation
%        5 = Full Truncation  (Least Bias)

%==============================
% Initialize Params/Vectors
%==============================

if nargin < 8
    scheme = 5;
end

Sigmav = SVModelParams.Sigmav;
v0     = SVModelParams.v0;
rho    = SVModelParams.rho;
theta = SVModelParams.theta;
eta   = SVModelParams.eta;


Spath = zeros(N_sim,M+1);
Spath(:,1) = S_0;

dt = T/M;
sqdt = sqrt(dt);
sqdtrho1 = sqdt*rho;
sqdtrho2 = sqdt*sqrt(1-rho^2);

Zeta = r - q;

edt = eta*dt;
coeffv = 1 - eta*dt;  %analytical drift for variance process
driftdt = eta*theta*dt;

vOld   = v0*ones(N_sim,1); %used to store variance process
sqvOld = sqrt(v0)*ones(N_sim,1);

if scheme == 1 %Absorbption
    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sqdt*sqvOld.*W1);  %log scheme

        vNew = driftdt + vOld*coeffv + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
        
        vOld = max(0,vNew); %least biased scheme
        sqvOld = sqrt(vOld);
    end
    
elseif scheme == 2   %Reflection

    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sqdt*sqvOld.*W1);  %log scheme

        vNew = driftdt + vOld*coeffv + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
        vOld = abs(vNew); %Always stays positive, the "reflection" scheme

        sqvOld = sqrt(vOld);
    end
    
elseif scheme == 3   % Higham and Mao

    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*sqvOld.^2)*dt + sqdt*sqvOld.*W1);  %log scheme

        vNew = driftdt + vOld*coeffv + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
        vOld = vNew; 
        
        sqvOld = sqrt(abs(vOld));
    end
    
elseif scheme == 4   % Partial Truncation

    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*sqvOld.^2)*dt + sqdt*sqvOld.*W1);  %log scheme

        vNew = driftdt + vOld*coeffv + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
        
        vOld = vNew; 
        sqvOld = sqrt(max(0, vOld));
    end
    
elseif scheme == 5   % Full Truncation

    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*sqvOld.^2)*dt + sqdt*sqvOld.*W1);  %log scheme

        vNew =  vOld  - edt*(max(0, vOld) - theta)  + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
        vOld = vNew;
        sqvOld =  sqrt(max(0, vOld)); 
    end

end

