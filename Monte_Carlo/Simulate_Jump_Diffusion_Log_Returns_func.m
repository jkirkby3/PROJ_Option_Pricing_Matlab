function Rets = Simulate_Jump_Diffusion_Log_Returns_func( N_sim, dt, r, q, sigma, jumpModel, jumpParams )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Simulates Log Returns of Jump Diffusion Models with jumps (including simple Black-Scholes, no jumps)
% Returns: paths of dimension (N_sim, M+1), since they include S_0 
%          ... Simulates N_sim paths, each row is a full path starting from S_0, ending with S_M (M+1 points in path)
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% N_sim = # paths
% M = #time steps on [0,T], time step is dt=T/M, so each path has M+1 points
% T = time to maturity, ie path is on [0,T]
% S_0 = initial underlying value (e.g. S_0=100)
% r = interst rate (e.g. r = 0.05)
% q = dividend yield (e.g. q = 0.05)
% sigma = diffusion parameter (e.g. sigma = 0.2)
%
%===================================
% jumpModel: 0 = NoJumps, 1 = NormalJumps, 2 = DEJumps, 3 = MixedNormalJumps
%===================================
% jumpParams = paramters container containing all necessary params for models,
%            : if jumpModel = 0, no jump params are needed
%            : if jumpModel > 0, jumpParams must contain lambda, kappa, and any other model specific params (see below)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

Sigsqdt = sigma*sqrt(dt);
drift = (Zeta - .5*sigma^2)*dt;
       
if jumpModel == 0  % Just a drifted Brownian motion
    Rets = drift + Sigsqdt*randn(N_sim,1); 
else
    Poi = PoissonRnd(N_sim, lamdt);  %Generate Poisson Column Vector of size N_Sim
    for n = 1:N_sim
        if Poi(n)>0
            Poi(n) = JumpFunc(Poi(n));  %JumpFunc(Poi(n)) sums up Poi(n) many jumps from the jump distribution
        end
    end

    Rets = drift + Poi + Sigsqdt*randn(N_sim,1); 
end


end

