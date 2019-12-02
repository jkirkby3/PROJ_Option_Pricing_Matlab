%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SWING OPTION PRICER (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Swing options in Levy Models using the PROJ method 
%               Supports 1) Constant and Linear Recovery Type Contracts
%                        2) Fixed Rights Contracts
%              
% Author:      Justin Kirkby
% References:  (1) Swing Option Pricing by Dynamic Programming with B-Spline Density Projection
%                   Int. J. Theoretical and App. Finance (2019)
%              (2) Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform, SIAM J. Financial Math., 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('./coeff_funcs')
addpath('../RN_CHF')
addpath('../Helper_Functions')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_0  = 40;  %Initial price

% Payoff constants (see paper for definition of payoff function)
K1  = 10;  %S_min
K2  = 20;  %K_d
K3  = 40;  %K_a
K4  = 50;  %S_max

r    = .05;  %Interest rate
q    = .00;  %dividend yield

% Common Contract Params
Dmax = 5;   % Max allowable swing size. Permissible swing size for any action is D={0,1,2,...,Dmax},
T    = 1;   % Time to maturity (in years)

% Contract Specific Params
swing_type = 3;  % 1 = Constant Recovery, 2 = Linear Recovery, 3 = Fixed Rights

if swing_type == 1
    tau1 = 1/4;   %%% Constant recovery time (in years) between swing actions tau_R(1)
    Mtau  = 12;  % Bermudan rights defintion (see paper)
elseif swing_type == 2
    rho_tau = 1/6;  % Linear Recovery time depends on size of exercise, this much of a year per unit (D) exercised
    Mtau  = 12;  % Bermudan rights defintion (see paper)
elseif swing_type == 3
    Ns = 5;  % Number of fixed rights
    M = 50;  % Number of Bermudan exercise dates available
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS  (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1;   %See Models Below (e.g. model 1 is Black Scholes), and choose specific params
params = {};

if model == 1 %BSM (Black Scholes Merton)
    params.sigmaBSM = 0.3;    %CHOOSE   
    
elseif model == 2 %CGMY
%     params.C  = 4; 
%     params.G  = 5; 
%     params.MM = 25; 
%     params.Y  = 0.8;

    params.C  = 1; 
    params.G  = 20; 
    params.MM = 30; 
    params.Y  = 1.7;

elseif model == 3 %NIG
    params.alpha = 15;
    params.beta  = -5;
    params.delta = 0.5;
    
elseif model == 4 %MJD (Merton Jump Diffusion)
    params.sigma  = 0.12;
    params.lam    = 0.4;
    params.muj    = -0.12;
    params.sigmaj = 0.18;
    
elseif model == 5 %Kou Double Expo
    params.sigma = 0.15;
    params.lam   = 3;
    params.p_up  = 0.2;
    params.eta1  = 25;
    params.eta2  = 10;
    
elseif model == 6 % Heston Model  
    params.v_0 = 0.0175; % initial variance
    params.theta = 0.0398;   % long term variance level
    params.kappa =1.5768;   % rate of variance mean reversion
    params.sigma_v = 0.5751;   % volatility of variance
    params.rho = -0.5711;   % correlation between Brownian motions
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 3) CHOOSE PROJ PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L1    = 20;  % sets gridwidth
logN  = 14;  % Uses N = 2^logN  gridpoints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelInput = getModelInput(model, T, r, q, params);
alpha = getTruncationAlpha(T, L1, modelInput, model);
N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]
Ks = [K1,K2,K3,K4];
T_0  = 0;  % Keep this at 0, this is inception. Minor modification needs to be made to price with T_0 > 0


if swing_type == 1  % Constant Recovery
    tic
    price = PROJ_Swing_ConstantRecovery_Aug( r,S_0,Dmax,T_0,T,tau1,Mtau,N,alpha,modelInput.rnSYMB,Ks);
    toc
elseif swing_type == 2  % Linear Recovery
    tic
    price = PROJ_Swing_LinearRec(r,S_0,Dmax,rho_tau,T_0,T,Mtau,N,alpha,modelInput.rnSYMB,Ks);
    toc
elseif swing_type == 3  % Fixed Rights
    modelInput = getModelInput(model, T/M, r, q, params);
    tic
    price = PROJ_Swing_FixedRights(M, N, alpha, modelInput.rnCHF, r, Dmax, T_0, T, Ns, S_0, Ks);
    toc
end

fprintf('%.8f \n', price)


