%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price European options under Stochastic Volatility Models (with Jumps) using Monte Carlo simulation
% Author:      Justin Kirkby
% For more details on these models see:
%          (1) A General Framework for discretely sampled realized
%              variance derivatives in stocahstic volatility models with
%              jumps, EJOR, 2017
%          (2) A unified approach to Bermudan and Barrier options under stochastic
%               volatility models with jumps. J. Economic Dynamics and Control, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);

addpath('../')

% ---------------------
%  Contract/Market Params
% ---------------------
call = 1;    %For call use 1 (else, its a put)
S_0  = 100;  %Initial price
r    = .05;  %Interest rate
q    = .00;  %dividend yield
T    = 1;    %Time (in years)
Kvec = S_0*[.85 .90 .95 1 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.5 1.6];   % strikes to price

%%%========================
%%%% Select Stochastic Volatility Model
%%%========================
model = 1;    % 1 = Heston (output compares with analytical)
              % 2 = Stein-Stein
              % 3 = 3/2 Model
              % 4 = 4/2 Model
              % 5 = Hull White (output compares with analytical)
              % 6 = Scott
              % 7 = Alpha-Hypergeometric

%%%========================
%%%% Select Jump Model
%%%========================
jumpModel = 0;    % 0 = No Jumps 
                  % 1 = Normal Jumps
                  % 2 = Double Exponential Jumps
                  % 3 = Mixed normal Jumps

% ---------------------
% Sim Params
% ---------------------
N_sim = 10^4;
M = 500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


jumpParams = {};

if jumpModel == 1 %Normal Jumps, e.g. Merton
    lambda = 1;  
    muJ = -.10; 
    sigJ = 0.3;
    
    jumpParams.kappa = exp(muJ + .5*sigJ^2)-1;  jumpParams.lambda = lambda; jumpParams.muJ = muJ; jumpParams.sigJ = sigJ;

elseif jumpModel == 2 %Double Exponenial Jumps     
    lambda = 1;
    p_up   = 0.5; % up jump probability    
    eta1   = 25;
    eta2   = 30;
    
    kappa  = p_up*eta1/(eta1-1)+(1-p_up)*eta2/(eta2+1)-1;
    jumpParams.lambda = lambda; jumpParams.kappa = kappa; jumpParams.eta1 = eta1; jumpParams.eta2 = eta2; jumpParams.p_up = p_up;    

elseif jumpModel == 3 %Mixed normal Jumps
    lambda = 1; 
    a1 = -0.05; 
    b1 = 0.07;    
    a2 = 0.02; 
    b2 = 0.03;
    p_up = 0.6;

    kappa = p_up*exp(a1 + .5*b1^2)+ (1-p_up)*exp(a2 + .5*b2^2)  -1;
    jumpParams.lambda = lambda; jumpParams.kappa = kappa; jumpParams.a1 = a1; jumpParams.b1 = b1; jumpParams.a2 = a2; jumpParams.b2 = b2; jumpParams.p_up = p_up;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Set the Stochastic Volatility Model Component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if model == 1
    %%%==============================
    %%% HESTON MODEL  Parameters
    %%%==============================
    modparam.eta    = 4;
    modparam.theta  = 0.035; 
    modparam.rho    = -0.75;
    modparam.Sigmav = 0.15;
    modparam.v0     = 0.04; 
    
elseif model == 2
    %%%=============================================================
    %%% STEIN-STEIN MODEL  Parameters
    %%%=============================================================
    modparam.eta    = 2; 
    modparam.theta  = 0.18; 
    modparam.Sigmav = 0.18; 
    modparam.v0     = 0.22; 
    modparam.rho    = -0.5; 
    
elseif model == 3
    %%%=============================================================
    %%% 3/2 MODEL  Parameters
    %%%=============================================================
    modparam.Sigmav = 0.10; 
    modparam.eta    = 3; 
    modparam.rho    = -0.7; 
    modparam.theta  = 0.04; 
    modparam.v0     = 0.04 ;
    
elseif model == 4
    %%%=============================================================
    %%% 4/2 MODEL  Parameters
    %%%=============================================================
    modparam.eta    = 3;
    modparam.theta  = 0.04; 
    modparam.rho    = -0.7; 
    modparam.Sigmav = 0.25; 
    modparam.v0     = 0.04; 
    modparam.aa     = 0.5; 
    modparam.bb     = 0.5*modparam.v0; 
    
elseif model == 5
    %%%=============================================================
    %%% HULL-WHITE MODEL  Parameters
    %%%=============================================================
    modparam.av     = 0.05; 
    modparam.rho    = -0.6;
    modparam.Sigmav = 0.6;
    modparam.v0     = 0.03; 
    
elseif model == 6
    %%%=============================================================
    %%% SCOTT MODEL  Parameters
    %%%============================================================= 
    modparam.eta    = 2; 
    modparam.theta  = log(0.16); 
    modparam.Sigmav = 0.20; 
    modparam.v0     = log(0.18);
    modparam.rho    = -0.9;

elseif model == 7
    %%%=============================================================
    %%% ALPHA-HYPERGEOMETRIC MODEL  Parameters
    %%%=============================================================
    modparam.rho    = -.9;
    modparam.Sigmav = .20; 
    modparam.v0     = log(0.17); 
    modparam.eta    = .05; 
    modparam.theta  = .2; 
    modparam.av     = 0.03;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Spath = Simulate_StochVol_Jumps_func( N_sim, M, T, S_0, r, q, model, modparam, jumpModel, jumpParams);
histogram(Spath(:,end))

disc = exp(-r*T);
[prices, stdErrs] = Price_MC_European_Strikes_func(Spath, disc, call, Kvec )

plot(Kvec, prices)
ylabel('price')
xlabel('strike')
grid on;


