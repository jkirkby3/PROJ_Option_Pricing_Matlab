%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Discrete Variance Swap / Option Pricer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Discrete Variance Swap / Options under stochastic volatility models (with jumps)
% Author:      Justin Kirkby
% References:  (1) A General Framework for discretely sampled realized
%              variance derivatives in stocahstic volatility models with
%              jumps, EJOR, 2017
%              (2) Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform, SIAM J. Financial Math., 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);

addpath('../Helper_Functions')
addpath('./Analytical_Swaps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contract = 1;     % Set contract = 1 for for variance swap, contract = 3 for variance call
K        = .00;   %strike, only matters for options, but still required for swap function
r        = 0.01;  % interest rate
q        = 0;     % dividend yield
T        = 0.5;     % Time to maturity
M        = 20;    % Number of monitoring dates on [0, T]

%%%========================
%%% PROJ Parameters
%%%========================
N  = 2^9;    % Number of basis elements
L1 = 14;     % Truncation width parameter
m_0 = 40;      % Number of variance states in CTMC
gamma = 5.5;     % State space grid width parameter  (5.5 is a decent starting point, but depends on the model)
varGridMult = .8;    % State space stretching paramter, chose in (0,1)... closer to 1 is more uniform grid
gridMethod = 4;     % Determines the grid method... 4 is a good choice

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Jump Model Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jumpParams = {};

if jumpModel == 0    %%%%% NO Jumps
    jumpParams.Nothing = 0;
    psi_J = @(u)0*[u>0];
    
    c2Jump = 0;
    c4Jump = 0;
    
elseif jumpModel == 1  %%%% Normal Jumps
    lambda = 1;  muJ = -.12;  sigJ = 0.15;
    
    jumpParams.kappa = exp(muJ + .5*sigJ^2)-1;
    jumpParams.lambda = lambda; jumpParams.muJ = muJ; jumpParams.sigJ = sigJ;
    psi_J = @(u) lambda*(exp(1i*u*muJ - .5*sigJ^2*u.^2)-1);    
    
    c2Jump = lambda*(muJ^2 +sigJ^2); %2nd cumulant of jump component
    c4Jump = lambda*(muJ^4 + 6*sigJ^2*muJ^2+3*sigJ^4*lambda);    
    
elseif jumpModel == 2 %%%% DE Jumps
    lambda = 1;
    p_up   = 0.5; % up jump probability    
    eta1   = 25;
    eta2   = 30;
    
    kappa  = p_up*eta1/(eta1-1)+(1-p_up)*eta2/(eta2+1)-1;
    jumpParams.lambda = lambda; jumpParams.kappa = kappa; jumpParams.eta1 = eta1; jumpParams.eta2 = eta2; jumpParams.p_up = p_up;     
    psi_J = @(u) lambda*(p_up*eta1./(eta1-1i*u) + (1-p_up)*eta2./(eta2+1i*u) -1) ;
    
    c2Jump = 2*lambda*p_up/eta1^2 + 2*lambda*(1-p_up)/eta2^2; %2nd cumulant of jump component
    c4Jump = 24*lambda*(p_up/eta1^4 + (1-p_up)/eta2^4);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Set the Stochastic Volatility Model Component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if model == 1
    %%%==============================
    %%% HESTON MODEL  Parameters
    %%%==============================
    modparam.eta    = 3.99;
    modparam.theta  = 0.014; 
    modparam.rho    = -0.79;
    modparam.Sigmav = 0.27;
    modparam.v0     = (.0994)^2; 
    
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PRICE CONTACT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alph = GetAlph_DisreteVariance( c2Jump, c4Jump, model, modparam, T, L1 );
PROJ_Price = PROJ_DiscreteVariance_StochVol( N,alph,M,r,T,K,m_0,psi_J, model, modparam, gridMethod, gamma, varGridMult, contract );
fprintf('PROJ Price: %.8f \n', PROJ_Price)

%%% In the special cases where analytic prices are known, also print the error
if model == 1 && jumpModel == 0 && contract == 1
   [ref, KcH] = hestonfairstrike(r, modparam.v0, modparam.theta, modparam.eta, modparam.Sigmav, T, modparam.rho, M);
   fprintf('Analytical Price: %.8f \n', ref)
   fprintf('Error: %.3e \n', PROJ_Price - ref)
   
elseif model == 5 && jumpModel == 0 && contract == 1
   [ref, KcH] = hullwhitefairstrike(r, modparam.v0, modparam.Sigmav, modparam.av, T, modparam.rho, M);
   Error1 = PROJ_Price - ref;
   fprintf('Analytical Price: %.8f \n', ref)
   fprintf('Error: %.3e \n', PROJ_Price - ref)
end
