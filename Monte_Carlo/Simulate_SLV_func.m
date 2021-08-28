function Spath = Simulate_SLV_func( N_sim, M, T, S_0, r, q, model, params)
% Simulates Paths of Stochastic Local Volatility Models (basic Euler Scheme for most cases) 
% N_sim = # paths
% M = #time steps on [0,T], ie dt =T/M   
% Note: returns paths of dimension (N_sim,M+1), since they include S_0
%
% SVLModel:                         (parameters)
%        1 = Heston:                alpha, v0, rho, theta, eta
%        2 = SABR:                  alpha, v0, rho, beta
%        3 = Shifted SABR:          alpha, v0, rho, beta, shift
%        4 = Quadratic SLV:         alpha, v0, rho, theta, eta, a, b, c
%        5 = TanHyp-Heston:         alpha, v0, rho, theta, eta, beta
%        6 = Heston-SABR:           alpha, v0, rho, theta, eta, beta
%
%==============================
% Initialize Common Params/Vectors
%==============================
alpha  = params.alpha;   % (vol of vol)
v0     = params.v0;  % initial vol/var
rho    = params.rho; % covar bewteen asset and vol innovations

Sigmav = alpha; %NOTE: alpha is same as Sigmav (vol of vol)

if model == 1 || model == 4 || model == 5 || model == 6
    theta = params.theta;
    eta   = params.eta;
end

Spath = zeros(N_sim,M+1);
Spath(:,1) = S_0;

dt = T/M;
sqdt = sqrt(dt);
sqdtrho1 = sqdt*rho;
sqdtrho2 = sqdt*sqrt(1-rho^2);

Zeta = r - q;
%==============================
% Simulate Based on Specific StochVol Model
%==============================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if model == 1 %HESTON MODEL
       
    expEta = exp(-eta*dt);
    driftv = theta*(1 - expEta);  %analytical drift for variance process
    vOld   = v0*ones(N_sim,1); %used to store variance process
    sqvOld = sqrt(v0)*ones(N_sim,1);
    
    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sqdt*sqvOld.*W1);  %log scheme
        %Spath(:,m+1) = Spath(:,m)*( (1+r)*dt + sqrt(vOld)*(sqrho1*W1 + sqrho2*W2));  %level scheme
        vNew = driftv + vOld*expEta + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
        %vNew = (sqvOld +.5*Sigmav*(sqdtrho1*W1 + sqdtrho2*W2)).^2 + eta*(theta - vOld)*dt -.25*Sigmav*dt; %Milstein
        vOld = abs(vNew); %Always stays positive, the "reflection" scheme
        %vOld = max(0,vNew); %least biased scheme
        sqvOld = sqrt(vOld);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif model == 2  % SABR 
    beta  = params.beta;
    cons1 = -.5*(alpha)^2*dt;
    vOld  = v0*ones(N_sim,1); %used to store variance process
    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        Spath(:,m+1) = max(0, Spath(:,m) + Spath(:,m).^beta.* vOld .*sqdt.*W1);  %level scheme
        vOld = vOld.*exp(cons1 + alpha*(sqdtrho1*W1 + sqdtrho2*W2));
        %vOld = vOld + alpha*vOld.*(sqdtrho1*W1 + sqdtrho2*W2);
    end

%     %%% LOG SCHEME 
%     cons2 = 2*(beta - 1);
%     for m = 1:M
%         W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
%         Spath(:,m+1) = Spath(:,m).*exp(-.5*vOld.^2.*Spath(:,m).^cons2*dt + Spath(:,m).^(beta-1).* vOld .*sqdt.*W1);  %level scheme
%         vOld = vOld.*exp(cons1 + alpha*(sqdtrho1*W1 + sqdtrho2*W2));
%     end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif model == 3  % SHIFTED SABR
    beta  = params.beta;
    shift = params.shift;
    cons1 = -.5*(alpha)^2*dt;
    vOld  = v0*ones(N_sim,1); %used to store variance process
    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        %%% Note: because of shift, -shift <= S_t
        Spath(:,m+1) = max(-shift, Spath(:,m) + (Spath(:,m) + shift).^beta.* vOld .*sqdt.*W1);  %level scheme
        vOld = vOld.*exp(cons1 + alpha*(sqdtrho1*W1 + sqdtrho2*W2));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif model == 4  % QUADRATIC SLV
    a = params.a; b = params.b; c = params.c;
    
    expEta = exp(-eta*dt);
    driftv = theta*(1 - expEta);  %analytical drift for variance process
    vOld   = v0*ones(N_sim,1); %used to store variance process
    sqvOld = sqrt(v0)*ones(N_sim,1);
    
    drift = 1+ Zeta*dt;
    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        
        Spath(:,m+1) = max(0,Spath(:,m)*drift + sqdt*sqvOld.*W1.*(Spath(:,m).*(a*Spath(:,m) +b) + c));
        
        %Spath(:,m+1) = Spath(:,m)*( (1+r)*dt + sqrt(vOld)*(sqrho1*W1 + sqrho2*W2));  %level scheme
        vNew = driftv + vOld*expEta + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
        %vNew = (sqvOld +.5*Sigmav*(sqdtrho1*W1 + sqdtrho2*W2)).^2 + eta*(theta - vOld)*dt -.25*Sigmav*dt; %Milstein
        %vOld = abs(vNew); %Always stays positive, the "reflection" scheme
        vOld = max(0,vNew); %least biased scheme
        sqvOld = sqrt(vOld);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif model == 5  % TAN-HYP-HESTON
    beta = params.beta;
    
    expEta = exp(-eta*dt);
    driftv = theta*(1 - expEta);  %analytical drift for variance process
    vOld   = v0*ones(N_sim,1); %used to store variance process
    sqvOld = sqrt(v0)*ones(N_sim,1);
    
    drift = 1+ Zeta*dt;
    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        
        Spath(:,m+1) = max(0,Spath(:,m)*drift + sqdt*sqvOld.*W1.*tanh(beta*Spath(:,m)));     
        vNew = driftv + vOld*expEta + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
        %vNew = (sqvOld +.5*Sigmav*(sqdtrho1*W1 + sqdtrho2*W2)).^2 + eta*(theta - vOld)*dt -.25*Sigmav*dt; %Milstein
        %vOld = abs(vNew); %Always stays positive, the "reflection" scheme
        vOld = max(0,vNew); %least biased scheme
        sqvOld = sqrt(vOld);
    end

elseif model == 6   % HESTON-SABR
    drift = 1+ Zeta*dt;
    expEta = exp(-eta*dt);
    driftv = theta*(1 - expEta);  %analytical drift for variance process
    vOld   = v0*ones(N_sim,1); %used to store variance process
    sqvOld = sqrt(v0)*ones(N_sim,1);
    
    beta  = params.beta;
    
    for m = 1:M
        W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
        
        %%%% LEVEL scheme
        Spath(:,m+1) = max(0, Spath(:,m)*drift + Spath(:,m).^beta.* sqvOld .*sqdt.*W1); 
        
%         %%%% LOG scheme
%         Spath(:,m+1) = Spath(:,m).*exp(-.5*vOld.^2.*Spath(:,m).^cons2*dt + Spath(:,m).^(beta-1).* sqvOld .*sqdt.*W1); 
        
        vNew = driftv + vOld*expEta + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
        vOld = abs(vNew); %Always stays positive, the "reflection" scheme
        %vOld = max(0,vNew); %least biased scheme
        sqvOld = sqrt(vOld);
    end


%     beta  = params.beta;
%     cons1 = -.5*(alpha)^2*dt;
%     vOld  = v0*ones(N_sim,1); %used to store variance process
%     for m = 1:M
%         W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
%         Spath(:,m+1) = max(0, Spath(:,m) + Spath(:,m).^beta.* vOld .*sqdt.*W1);  %level scheme
%         vOld = vOld.*exp(cons1 + alpha*(sqdtrho1*W1 + sqdtrho2*W2));
%         %vOld = vOld + alpha*vOld.*(sqdtrho1*W1 + sqdtrho2*W2);
%     end
%     
%     expEta = exp(-eta*dt);
%     driftv = theta*(1 - expEta);  %analytical drift for variance process
%     vOld   = v0*ones(N_sim,1); %used to store variance process
%     sqvOld = sqrt(v0)*ones(N_sim,1);
%     
%     for m = 1:M
%         W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
%         Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sqdt*sqvOld.*W1);  %log scheme
%         %Spath(:,m+1) = Spath(:,m)*( (1+r)*dt + sqrt(vOld)*(sqrho1*W1 + sqrho2*W2));  %level scheme
%         vNew = driftv + vOld*expEta + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
%         %vNew = (sqvOld +.5*Sigmav*(sqdtrho1*W1 + sqdtrho2*W2)).^2 + eta*(theta - vOld)*dt -.25*Sigmav*dt; %Milstein
%         vOld = abs(vNew); %Always stays positive, the "reflection" scheme
%         %vOld = max(0,vNew); %least biased scheme
%         sqvOld = sqrt(vOld);
%     end
   
end







% %%%% REFORMULATED VERSION
% etahat = eta*theta; thetahat = (eta + Sigmav^2)/etahat; Sigmavhat = -Sigmav; v0hat = 1/v0;
% vOld = v0hat*ones(N_sim,1); %used to store variance process
% sqvOld = sqrt(v0hat)*ones(N_sim,1);
% 
% expkt = exp(-etahat*dt);
% driftv = thetahat*(1 - expkt);  %analytical drift for variance process
% sqtSv = sqdt*Sigmavhat;  %analytical drift for variance process
% 
% for m = 1:M
%     W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
% 
%     Spath(:,m+1) = Spath(:,m).*exp((r - .5*vOld)*dt + (sqrho1*W1 + sqrho2*W2)./sqvOld);  %log scheme
%     %Spath(:,m+1) = Spath(:,m)*( (1+r)*dt + sqrt(vOld)*(sqrho1*W1 + sqrho2*W2));  %level scheme
% 
%     vNew = driftv + vOld*expkt + sqtSv*sqvOld.*W1;
% 
%     vOld = abs(vNew); %Always stays positive, the "reflection" scheme
%     sqvOld = sqrt(vOld);
% end
