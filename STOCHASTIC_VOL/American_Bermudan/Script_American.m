%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Script For American Put Options under Stochastic Volatility,
% Using Regime Switching Fourier Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../Helper_Functions')

%%%-------------------------------
S_0 = 50;
W   = 55;  %strike
r   = 0.05; 
T   = .25;
M   = 80;

%%%----------------------------
N    = 2^10;    %number of points in density expansion... Value grid size is K:=N/2
alph = 6;  %density projection grid on [-alpha,alpha]
%%%----------------------------
m_0           = 30;
gamma         = 3.3;
gridMethod    = 4;
gridMultParam = 0.2;

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


if model == 1 
    %%%=============================================================
    %%% HESTON  Parameters
    %%%=============================================================
    %etaH =4; thetaH =0.035; rhoH = -0.6; SigmavH =0.15; v0H=0.04;  %%Test Set 1
    %etaH =1.5; thetaH =0.035; rhoH = -0.8; SigmavH =0.12; v0H=0.04; %%Test Set 2
    %etaH =3; thetaH =0.04; rhoH = -0.1; SigmavH =0.1; v0H=0.04;  %%Test Set 3

    etaH =4; thetaH =0.035; rhoH = -0.75; SigmavH =0.15; v0H=0.04;

    %%%%%%%%%%%%%%%%%%%%%%%%%    
    modparam.eta = etaH; modparam.theta = thetaH; modparam.rho = rhoH; modparam.Sigmav = SigmavH; modparam.v0 = v0H;
    %%%%%%%%%%%%%%%%%%%%%%%%%

elseif model == 2 
    %%%=============================================================
    %%% STEIN-STEIN (ie Schoebel-Zhu)  Parameters
    %%%=============================================================     
    etaS = 3; thetaS = 0.18; SigmavS = 0.15; v0S = 0.15; rhoS = -0.6;  %%Test Set 1
    %etaS = 2; thetaS = 0.18; SigmavS = 0.18; v0S = 0.22; rhoS = -0.5; %%Test Set 2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    modparam.eta = etaS; modparam.theta = thetaS; modparam.rho = rhoS; modparam.Sigmav = SigmavS; modparam.v0 = v0S;
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    
elseif model == 3 
    %%=============================================================
    %%% 3/2 MODEL  Parameters
    %%=============================================================
    Sigmav32 =0.15; eta32 = 4; rho32= -0.6; theta32 =0.03;  v032 =0.03 ;  %Test Set 1
    %Sigmav32 =0.12; eta32 = 2; rho32= -.3; theta32 =0.04;  v032 = 0.035;  %Test Set 2
    %Sigmav32 =0.10; eta32 = 3; rho32= -0.7; theta32 =0.04;  v032 =0.04 ;%Test Set 3
    
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    modparam.eta = eta32; modparam.theta = theta32; modparam.rho = rho32; modparam.Sigmav = Sigmav32; modparam.v0 = v032;
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    
elseif model == 4 
    %%%=============================================================
    %%% 4/2 MODEL  Parameters
    %%%=============================================================
    eta42 = 0.5; theta42 = 0.035; rho42=-0.6; Sigmav42=0.15; v042 = 0.04; aa = 0.5; bb = 0.5*v042;  %%Test Set 1
    %eta42 = 1.8; theta42 = 0.04; rho42=-0.7; Sigmav42=0.1; v042 = 0.04; aa = 0.5; bb = 0.5*v042;  %%Test Set 2
    %eta42 = 3; theta42 = 0.04; rho42=-0.1; Sigmav42=0.1; v042 = 0.04; aa = 0.5; bb = 0.5*v042;  %%Test Set 3

    %%%%%%%%%%%%%%%%%%%%%%%%%  
    modparam.eta = eta42; modparam.theta = theta42; modparam.rho = rho42; modparam.Sigmav = Sigmav42; modparam.v0 = v042;
    modparam.aa = aa; modparam.bb = bb;
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    
elseif model == 5 
    %%%=============================================================
    %%% HULL-WHITE  Parameters
    %%%=============================================================
    av = 0.01; rhoHW = -0.6; SigmavHW = 0.25; v0HW = 0.04;   %%Test Set 1
    %av = 0.03; rhoHW = -0.7; SigmavHW = 0.10; v0HW = 0.03;  %%Test Set 2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    modparam.rho = rhoHW; modparam.Sigmav = SigmavHW; modparam.v0 = v0HW; modparam.av = av;
    %%%%%%%%%%%%%%%%%%%%%%%%%     

elseif model == 6 
    % %%%=============================================================
    % %%% SCOTT  Parameters
    % %%%=============================================================   
    etaSc = 3; thetaSc = log(0.15); SigmavSc = 0.15; v0Sc = log(0.15); rhoSc = -.6;  %%Test Set 1
    %etaSc = 2; thetaSc = log(0.16); SigmavSc = 0.20; v0Sc = log(0.18); rhoSc = -.9;  %%Test Set 2
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
price = American_StochVol_func( N,alph,M,r,T,S_0,W,m_0,psi_J,model, modparam, gridMethod, gamma, gridMultParam);
toc
fprintf('%.8f \n', price)


