[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);

addpath('../Helper_Functions')
addpath('./Analytical_Swaps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to generate prices for variance swaps and options in SV models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stochastic Volatility Model Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if model ==1
    %%%=============================================================
    %%% HESTON MODEL  Parameters
    %%%=============================================================

    %etaH =3.99; thetaH =0.014; rhoH = -0.79; SigmavH =0.27; v0H=(.0994)^2;
    etaH =3.99; thetaH =0.014; rhoH = -0.79; SigmavH =0.27; v0H=(.0994)^2;
%     etaH = 3;
%     rhoH = 0;
%     thetaH = 0.04;
%     SigmavH = 0.1;
%     v0H = 0.04;
    
    % % % etaH =6.21; thetaH =0.019; rhoH = -0.70; SigmavH =0.31; v0H=(.1011)^2;
    % % % etaH =3; thetaH =0.04; rhoH = -0.7; SigmavH =0.25; v0H=.03;  
    % % % etaH =1; thetaH =0.02; rhoH = -0.9; SigmavH =0.15; v0H=.04;  
    % % % etaH =0.5; thetaH =0.03; rhoH = -0.8; SigmavH =0.09; v0H=.04;  
    % % % etaH = 6.52; thetaH =0.0352; rhoH = -0.771; SigmavH =0.4601; v0H=.03;  
    % % % etaH =6.21; thetaH =0.019; rhoH = -0.7; SigmavH =0.31; v0H=(.1011)^2;     
    % % % etaH =3; thetaH =0.04; rhoH = -0.1; SigmavH =0.25; v0H=0.03;  
    % % % etaH = 3; thetaH =0.04; rhoH = -0.7; SigmavH =0.25; v0H=.04;
    % % % etaH = 6.5; thetaH =0.035; rhoH = -0.7; SigmavH =0.45; v0H=.03;  
    % % % etaH = 2.5; thetaH =0.03; rhoH = -0.7; SigmavH =0.20; v0H=.04; 
    % % % etaH = 1.5; thetaH =0.035; rhoH = -0.8; SigmavH =0.25; v0H=.045;  
    % % % etaH = 3; thetaH =0.045; rhoH = -0.6; SigmavH =0.3; v0H=.03;  

    %%%%%%%%%%%%
    modparam.eta = etaH; modparam.theta = thetaH; modparam.rho = rhoH; modparam.Sigmav = SigmavH; modparam.v0 = v0H;
    %%%%%%%%%%%

elseif model == 2
    %%%=============================================================
    %%% STEIN-STEIN MODEL  Parameters
    %%%=============================================================
    %etaS = 3; thetaS = 0.18; SigmavS = 0.15; v0S = 0.15; rhoS = -0.6;  %%Test Set 1
    etaS = 2; thetaS = 0.18; SigmavS = 0.18; v0S = 0.22; rhoS = -0.5; %%Test Set 2
    %etaS = 5; thetaS = 0.18; SigmavS = 0.15; v0S = 0.15; rhoS = -0.7; %%Test Set 3
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    modparam.eta = etaS; modparam.theta = thetaS; modparam.rho = rhoS; modparam.Sigmav = SigmavS; modparam.v0 = v0S;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif model == 3
    %%%=============================================================
    %%% 3/2 MODEL  Parameters
    %%%=============================================================
    % Sigmav32 =0.15; eta32 = 4; rho32= -0.6; theta32 =0.03;  v032 =0.03 ;  
    % Sigmav32 =0.12; eta32 = 2; rho32= -.3; theta32 =0.04;  v032 = 0.035; 
    Sigmav32 =0.10; eta32 = 3; rho32= -0.7; theta32 =0.04;  v032 =0.04 ; 

    %%%%%%%%%%%%%%%%%%%%%%%%%
    modparam.eta = eta32; modparam.theta = theta32; modparam.rho = rho32; modparam.Sigmav = Sigmav32; modparam.v0 = v032;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif model == 4  
    %%%=============================================================
    %%% 4/2 MODEL  Parameters
    %%%=============================================================
    eta42 = 3; theta42 = 0.04; rho42= -0.7; Sigmav42=0.25; v042 = 0.04; aa = 0.5; bb = 0.5*v042;
    %eta42 = 4; theta42 = 0.035; rho42=-0.7; Sigmav42=0.20; v042 = 0.04; aa =0.5; bb = 0.5*v042;  
    %eta42 = 3; theta42 = 0.04; rho42=-0.1; Sigmav42=0.15; v042 = 0.045; aa =0.5; bb = 0.5*v042; 
    %eta42 = 1.8; theta42 = 0.04; rho42=-0.7; Sigmav42=0.2; v042 = 0.04; aa =0.5; bb = 0.5*v042; 

    %%%%%%%%%%%%%%%%%%%%%%%%%
    modparam.eta = eta42; modparam.theta = theta42; modparam.rho = rho42; modparam.Sigmav = Sigmav42; modparam.v0 = v042;
    modparam.aa = aa; modparam.bb = bb;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif model == 5  
    %%%=============================================================
    %%% HULL-WHITE MODEL  Parameters
    %%%=============================================================
    %av = 0.05; rhoHW = -0.7; SigmavHW = 0.6; v0HW = 0.05; r = 0.03; %%Test Set 3
    %av = 0.01; rhoHW = -0.6; SigmavHW = 0.25; v0HW = 0.03; r = 0.03; %%Test Set 4
    %av = 0.001; rhoHW = -0.9; SigmavHW = 0.35; v0HW = 0.03; r = 0.03;
    av = 0.09; rhoHW = -0.2; SigmavHW = 0.4; v0HW = 0.034; r = 0.05;

    %%%%%%%%%%%%%%%%%%%%%%%%%   %NOTE: we used "mu" and "av" interchangibly for Hull white
    modparam.rho = rhoHW; modparam.Sigmav = SigmavHW; modparam.v0 = v0HW; modparam.av = av;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif model == 6  
    %%%=============================================================
    %%% SCOTT MODEL  Parameters
    %%%============================================================= 
    %etaSc = 3; thetaSc = log(0.15); SigmavSc = 0.15; v0Sc = log(0.15); rhoSc = -.6;  %%Test Set 1
    etaSc = 2; thetaSc = log(0.16); SigmavSc = 0.20; v0Sc = log(0.18); rhoSc = -.9;  %%Test Set 2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    modparam.eta = etaSc; modparam.theta = thetaSc; modparam.rho = rhoSc; modparam.Sigmav = SigmavSc; modparam.v0 = v0Sc;
    %%%%%%%%%%%%%%%%%%%%%%%%%

elseif model == 7 
    %%%=============================================================
    %%% ALPHA-HYPER MODEL  Parameters
    %%%=============================================================
    rhoAH = -.6; SigmavAH = .25; v0AH = log(0.19); etaAH = .01; thetaAH = .04; avAH = 0.03; %Test Set 1
    %rhoAH = -.9; SigmavAH = .20; v0AH = log(0.17); etaAH = .05; thetaAH = .2; avAH = 0.03; % Test Set 2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    modparam.rho = rhoAH;  modparam.Sigmav = SigmavAH;  modparam.v0 = v0AH; modparam.eta = etaAH; modparam.av = avAH; modparam.theta = thetaAH;
    %%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PRICE CONTACT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alph = GetAlph_DisreteVariance( c2Jump, c4Jump, model, modparam, T, L1 );
PROJ_Price = DiscreteVariance_PROJ( N,alph,M,r,T,K,m_0,psi_J, model, modparam, gridMethod, gamma, varGridMult, contract );
fprintf('PROJ Price: %.8f \n', PROJ_Price)

%%% In the special cases where analytic prices are known, also print the error
if model == 1 && jumpModel == 0 && contract == 1
   [ref, KcH] = hestonfairstrike(r, v0H, thetaH, etaH, SigmavH, T, rhoH, M);
   fprintf('Analytical Price: %.8f \n', ref)
   fprintf('Error: %.3e \n', PROJ_Price - ref)
   
elseif model == 5 && jumpModel == 0 && contract == 1
   [ref, KcH] = hullwhitefairstrike(r, v0HW, SigmavHW, av, T, rhoHW, M);
   Error1 = PROJ_Price - ref;
   fprintf('Analytical Price: %.8f \n', ref)
   fprintf('Error: %.3e \n', PROJ_Price - ref)
end
