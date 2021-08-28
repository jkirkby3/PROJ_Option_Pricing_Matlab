function Spath = Simulate_StochVol_Jumps_func( N_sim, M, T, S_0, r, q, SVModel, SVModelParams, jumpModel, jumpParams)
% Simulates Paths of Stochastic Volatility Models (basic Euler Scheme for most cases) with jumps
% By default, if no jumpModel or jumpParams are passed, it defaults to stochastic vol without jumps
%
% N_sim = # paths
% M = #time steps on [0,T], ie dt =T/M   
% T = final time
% S_0 = initial underlying (spot) value
% r = interest rate
% q = dividend / convenience yield
%
% Note: returns paths of dimension (N_sim,M+1), since they include S_0
%
%===================================
% jumpModel: 0 = NoJumps, 1 = NormalJumps, 2 = DEJumps, 3 = MixedNormalJumps
%===================================
% SVModel:    (with parameters)
%        1 = HESTON:       Sigmav, v0, rho, eta, theta
%        2 = STEIN-STEIN:  Sigmav, v0, rho, eta, theta
%        3 = 3/2 MODEL:    Sigmav, v0, rho, eta, theta
%        4 = 4/2 MODEL:    Sigmav, v0, rho, eta, theta, aa, bb
%        5 = HULL-WHITE:   Sigmav, v0, rho, av
%        6 = SCOTT:        Sigmav, v0, rho, eta, theta
%        7 = ALPHA-HYPER:  Sigmav, v0, rho, eta, theta
%        8 = "VAR" MODEL:  Sigmav, v0, rho, eta, theta
%        9 = Jacobi Model: vmin, vmax, Sigmav, v0, rho, eta, theta
%
%==============================
% Initialize Common Params/Vectors
%==============================
Sigmav = SVModelParams.Sigmav;
v0     = SVModelParams.v0;
rho    = SVModelParams.rho;

if nargin < 10  % Defaults to the case of No Jumps
    jumpModel = 0; jumpParams={};
end

if SVModel == 1 || SVModel == 2 || SVModel == 3 || SVModel == 4 || SVModel == 6 || SVModel == 7 || SVModel == 8 || SVModel == 9
    theta = SVModelParams.theta;
    eta   = SVModelParams.eta;
end
if SVModel == 4 %4/2
    aa = SVModelParams.aa;
    bb = SVModelParams.bb;
end
if SVModel == 5 || SVModel == 7 %hull-white / alpha-hyper
    av = SVModelParams.av;
end

Spath = zeros(N_sim,M+1);
Spath(:,1) = S_0;

dt = T/M;
sqdt = sqrt(dt);
sqdtrho1 = sqdt*rho;
sqdtrho2 = sqdt*sqrt(1-rho^2);

%==============================
% Initialize Jump Model Params and JumpFunc (function handle)
%==============================
%%% NOTE:  Jump Model is of the form in LOG space
%%% X(m+1) = X(m) + drift + Brownian Component + sum(Jumps on [m,m+1])
%%% By Jump we mean log(Y), e.g. in Merton Model, Jump ~ Normal (since we are in log space )

if jumpModel > 0 %ie if there are jumps in the model
    lambda = jumpParams.lambda;
    kappa  = jumpParams.kappa;

    Zeta = r - q - lambda*kappa;  %NOTE: we are redefining r to include compensation for jump component
    lamdt = lambda*dt;
    
    if jumpModel == 1 %Normal Jumps, e.g. Merton
        muJ  = jumpParams.muJ;
        sigJ = jumpParams.sigJ;
        JumpFunc = @(n) sum(muJ +sigJ*randn(n,1)); %Generates n independent jumps and sums them
        
    elseif jumpModel == 2 %Double Exponenial Jumps     
        p_up = jumpParams.p_up;       
        eta1 = jumpParams.eta1;
        eta2 = jumpParams.eta2;       
        JumpFunc = @(n) sum(DoubleExpoRnd(n,p_up, eta1,eta2));
        
    elseif jumpModel == 3 %Mixed normal Jumps
        p_up = jumpParams.p_up;
        a1 = jumpParams.a1;  b1 = jumpParams.b1;
        a2 = jumpParams.a2;  b2 = jumpParams.b2;
        JumpFunc = @(n) sum(MixedNormalRnd(n, p_up, a1,b1,a2,b2));
    end
else 
    Zeta = r - q;
end
%==============================
% Simulate Based on Specific StochVol Model
%==============================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if SVModel == 1 %HESTON MODEL
    % NOTE: Uses Full Truncation Scheme studied in Lord et. al. (2010)
       
    vOld   = v0*ones(N_sim,1); %used to store variance process
    sqvOld = sqrt(v0)*ones(N_sim,1);
    edt = eta*dt;
    
    if jumpModel == 0
        for m = 1:M
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*sqvOld.^2)*dt + sqdt*sqvOld.*W1);  %log scheme

            vNew =  vOld  - edt*(max(0, vOld) - theta)  + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
            vOld = vNew;
            sqvOld =  sqrt(max(0, vOld)); 
            
        end
    else
        for m = 1:M
            Poi = PoissonRnd(N_sim, lamdt);  %Generate Poisson Column Vector of size N_Sim
            sumJumpsVec = zeros(N_sim,1);
            for n = 1:N_sim
                if Poi(n)>0
                    sumJumpsVec(n) = JumpFunc(Poi(n));  %JumpFunc(Poi(n)) sums up Poi(n) many jumps from the jump distribution
                end
            end

            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*sqvOld.^2)*dt + sumJumpsVec + sqdt*sqvOld.*W1);  %log scheme

            vNew =  vOld  - edt*(max(0, vOld) - theta)  + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
            vOld = vNew;
            sqvOld =  sqrt(max(0, vOld)); 
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif SVModel == 2  % STEIN-STEIN
    driftv = eta*theta*dt;
    consv = 1-eta*dt;
    %conss = 1+r*dt;
    vOld = v0*ones(N_sim,1); %used to store variance process
    
    if jumpModel == 0
        for m = 1:M
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld.^2)*dt + sqdt*vOld.*W1);  %log scheme
            %Spath(:,m+1) = Spath(:,m).*(conss + vOld.*(sqrho1*W1 + sqrho2*W2));  %level scheme
            vNew = driftv + vOld*consv + Sigmav*(sqdtrho1*W1 + sqdtrho2*W2);
            %vOld = abs(vNew); %Always stays positive, the "reflection" scheme
            vOld = max(0,vNew); %least biased scheme
        end
    else
        for m = 1:M
            Poi = PoissonRnd(N_sim, lamdt);  %Generate Poisson Column Vector of size N_Sim
            %Poi = poissrnd(lamdt, N_sim,1);  %Generate Poisson Column Vector of size N_Sim
            sumJumpsVec = zeros(N_sim,1);
            for n = 1:N_sim
                if Poi(n)>0
                    sumJumpsVec(n) = JumpFunc(Poi(n));  %JumpFunc(Poi(n)) sums up Poi(n) many jumps from the jump distribution
                end
            end
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld.^2)*dt + sumJumpsVec + sqdt*vOld.*W1);  %log scheme
            vNew = driftv + vOld*consv + Sigmav*(sqdtrho1*W1 + sqdtrho2*W2);
            %vOld = abs(vNew); %Always stays positive, the "reflection" scheme
            vOld = max(0,vNew); %least biased scheme
        end
    end 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif SVModel == 3 % 3/2 MODEL
      
    vOld   = v0*ones(N_sim,1); %used to store variance process
    sqvOld = sqrt(v0)*ones(N_sim,1);
    
    if jumpModel == 0  %%% 3/2 NO JUMPS
        for m = 1:M
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sqdt*sqvOld.*W1);  %log scheme
            %Spath(:,m+1) = Spath(:,m)*( (1+r)*dt + sqrt(vOld)*(sqrho1*W1 + sqrho2*W2));  %level scheme
            vNew = vOld.*( 1 + eta*(theta - vOld)*dt + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2));
            %vOld = abs(vNew); %Always stays positive, the "reflection" scheme
            vOld = max(0,vNew); %least biased scheme
            sqvOld = sqrt(vOld);
        end
    else %%% 3/2 WITH JUMPS
        for m = 1:M
            Poi = PoissonRnd(N_sim, lamdt);  %Generate Poisson Column Vector of size N_Sim
            %Poi = poissrnd(lamdt, N_sim,1);  %Generate Poisson Column Vector of size N_Sim
            sumJumpsVec = zeros(N_sim,1);
            for n = 1:N_sim
                if Poi(n)>0
                    sumJumpsVec(n) = JumpFunc(Poi(n));  %JumpFunc(Poi(n)) sums up Poi(n) many jumps from the jump distribution
                end
            end
            
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sumJumpsVec +  sqdt*sqvOld.*W1);  %log scheme
            vNew = vOld.*( 1 + eta*(theta - vOld)*dt + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2));
            %vOld = abs(vNew); %Always stays positive, the "reflection" scheme
            vOld = max(0,vNew); %least biased scheme
            sqvOld = sqrt(vOld);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
elseif SVModel == 4 % 4/2 MODEL
    expEta = exp(-eta*dt);
    driftv = theta*(1 - expEta);  %analytical drift for variance process

    vOld = v0*ones(N_sim,1); %used to store variance process
    sqvOld = sqrt(v0)*ones(N_sim,1);
    if jumpModel == 0
        for m = 1:M
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sqdt*(aa*sqvOld + bb./sqvOld).*W1);  %log scheme
            %Spath(:,m+1) = Spath(:,m)*( (1+r)*dt + sqrt(vOld)*(sqrho1*W1 + sqrho2*W2));  %level scheme
            vNew = driftv + vOld*expEta + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
            %vNew = (sqvOld +.5*Sigmav*sqdt*W1).^2 + eta*(theta - vOld)*dt -.25*Sigmav*dt; %Milstein
            vOld = abs(vNew); %Always stays positive, the "reflection" scheme       
            %vOld = max(0,vNew);
            sqvOld = sqrt(vOld);
        end
    else
        for m = 1:M
            Poi = PoissonRnd(N_sim, lamdt);  %Generate Poisson Column Vector of size N_Sim
            sumJumpsVec = zeros(N_sim,1);
            for n = 1:N_sim
                if Poi(n)>0
                    sumJumpsVec(n) = JumpFunc(Poi(n));  %JumpFunc(Poi(n)) sums up Poi(n) many jumps from the jump distribution
                end
            end          
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sumJumpsVec + sqdt*(aa*sqvOld + bb./sqvOld).*W1);  %log scheme
            vNew = driftv + vOld*expEta + Sigmav*sqvOld.*(sqdtrho1*W1 + sqdtrho2*W2);
            vOld = abs(vNew); %Always stays positive, the "reflection" scheme
            %vOld = max(0,vNew); %least biased scheme
            sqvOld = sqrt(vOld);
        end
    end 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

elseif SVModel == 5 % HULL-WHITE   

    driftv = (av - Sigmav^2/2)*dt;  %analytical drift for variance process
    vOld   = v0*ones(N_sim,1); %used to store variance process

    if jumpModel == 0
        for m = 1:M
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sqdt*sqrt(vOld).*W1);  %log scheme
            %Spath(:,m+1) = Spath(:,m)*( (1+r)*dt + sqrt(vOld)*(sqrho1*W1 + sqrho2*W2));  %level scheme
            vNew = vOld.*exp(driftv + Sigmav*(sqdtrho1*W1 + sqdtrho2*W2));    
            vOld = max(0,vNew); %least biased scheme
        end
    else
        for m = 1:M
            Poi = PoissonRnd(N_sim, lamdt);  %Generate Poisson Column Vector of size N_Sim
            sumJumpsVec = zeros(N_sim,1);
            for n = 1:N_sim
                if Poi(n)>0
                    sumJumpsVec(n) = JumpFunc(Poi(n));  %JumpFunc(Poi(n)) sums up Poi(n) many jumps from the jump distribution
                end
            end          
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sumJumpsVec + sqdt*sqrt(vOld).*W1);  %log scheme
            vNew = vOld.*exp(driftv + Sigmav*(sqdtrho1*W1 + sqdtrho2*W2));    
           vOld = max(0,vNew); %least biased scheme
        end
    end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif SVModel == 6 % SCOTT
    driftv = eta*theta*dt;
    consv = 1-eta*dt;
    vOld = v0*ones(N_sim,1); %used to store variance process
    
    if jumpModel == 0
        for m = 1:M
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*exp(2*vOld))*dt + sqdt*exp(vOld).*W1);  %log scheme
            vNew = driftv + vOld*consv + Sigmav*(sqdtrho1*W1 + sqdtrho2*W2);
            vOld = vNew; %Scott allowed to be negative    
        end
    else
        for m = 1:M
            Poi = PoissonRnd(N_sim, lamdt);  %Generate Poisson Column Vector of size N_Sim
            sumJumpsVec = zeros(N_sim,1);
            for n = 1:N_sim
                if Poi(n)>0
                    sumJumpsVec(n) = JumpFunc(Poi(n));  %JumpFunc(Poi(n)) sums up Poi(n) many jumps from the jump distribution
                end
            end          
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*exp(2*vOld))*dt +sumJumpsVec+ sqdt*exp(vOld).*W1);  %log scheme
            vNew = driftv + vOld*consv + Sigmav*(sqdtrho1*W1 + sqdtrho2*W2);
            vOld = vNew; %Scott allowed to be negative
        end
    end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
elseif SVModel == 7 % ALPHA-HPYERGEOMETRIC

    vOld = v0*ones(N_sim,1); %used to store variance process
    expvOld = exp(vOld);
    if jumpModel == 0
        for m = 1:M
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*expvOld.^2)*dt + sqdt*expvOld.*W1);  %log scheme
            vNew = vOld +(eta -theta*expvOld.^av)*dt + Sigmav*(sqdtrho1*W1 + sqdtrho2*W2);
            vOld = vNew; %alpha-h allowed to be negative
            expvOld = exp(vOld);
        end
    else
        for m = 1:M
            Poi = PoissonRnd(N_sim, lamdt);  %Generate Poisson Column Vector of size N_Sim
            sumJumpsVec = zeros(N_sim,1);
            for n = 1:N_sim
                if Poi(n)>0
                    sumJumpsVec(n) = JumpFunc(Poi(n));  %JumpFunc(Poi(n)) sums up Poi(n) many jumps from the jump distribution
                end
            end          
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*expvOld.^2)*dt + sumJumpsVec + sqdt*expvOld.*W1);  %log scheme
            vNew = vOld +(eta -theta*expvOld.^av)*dt + Sigmav*(sqdtrho1*W1 + sqdtrho2*W2);
            vOld = vNew; %alpha-h allowed to be negative
            expvOld = exp(vOld);
        end
    end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
elseif SVModel == 8 % VAR Model: See christoferrsen, Jacobs, Mimouni 2010
    if jumpModel == 0

    else
            %%%% ADD CODE FOR JUMP MODELS!!!!!!!!!!!!!!!!!!!
    end   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif SVModel == 9 % JACOBI Model
    vmin = SVModelParams.vmin;
    vmax = SVModelParams.vmax;
    denom = (sqrt(vmax) - sqrt(vmin))^2;
    driftv = eta*theta*dt;
    cons5 = (1-eta*dt);
    rhoSq = rho^2;
    vOld = v0*ones(N_sim,1); %used to store variance process
    QvOld = (vOld - vmin).*(vmax - vOld)/denom;
    if jumpModel == 0
        for m = 1:M
            W1 = randn(N_sim,1); W2 = randn(N_sim,1);  %Generate two Brownian motions
            Spath(:,m+1) = Spath(:,m).*exp((Zeta - .5*vOld)*dt + sqdt*sqrt(vOld - rhoSq*QvOld).*W1 + sqdtrho1*sqrt(QvOld).*W2 );  %log scheme

            vNew = driftv + cons5*vOld + Sigmav*sqdt*sqrt(QvOld).*W2 ;
            %vNew = vOld + eta*(theta - vOld)*dt + Sigmav*sqdt*sqrt(QvOld).*W2 ;
            
            %%% Variance process must lie between [vmin, vmax]
            vOld = min(vmax, max(vmin,vNew));  
            QvOld = (vOld - vmin).*(vmax - vOld)/denom;
            
        end        
    else
            %%%% ADD CODE FOR JUMP MODELS!!!!!!!!!!!!!!!!!!!
    end   
    
end



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
