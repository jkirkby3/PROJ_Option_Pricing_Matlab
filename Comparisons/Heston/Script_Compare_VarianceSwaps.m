%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VARIANCE SWAP COMPARISON (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Compare Methods For Variance Swaps Under Heston's Model
% Author:      Justin Kirkby
% 
% Methods: 1) Analytical 
%            - (Bernard and Cui 2014, Prices and Asymptotics for discrete variance swaps) 
%          2) Time-Changed Markov Chain 
%             - (Cui, Kirkby, Nguyen, 2019, A general framework for time-changed Markov processes and applications)
%             - Assumes rho = 0
%          3) SV-CTMC approximation 
%             - (A General Framework for discretely sampled realized
%              variance derivatives in stocahstic volatility models with
%              jumps, EJOR, Cui, Kirkby, Nguyen, 2017)
%          4) Monte Carlo
%               ... More to come
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r    = .02;  %Interest rate (NOTE: set to zero for comparison with Kahl-Jackel-Lord, based on Forward price)
q    = .00;  %dividend yield
T    = 1;    %Time (in years)
M    = 252;  %number of discrete monitoring points


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS  (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta = 2;  
theta = 0.03;
Sigmav = 0.3;
v0 = 0.04;
rho = 0.0;  % NOTE: when rho is not zero, Time-Changed CTMC becomes approximation using rho=0 (which is often very good)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = {};

params.v0 = v0; % initial variance
params.theta = theta;   % long term variance level
params.eta = eta;   % rate of variance mean reversion
params.Sigmav = Sigmav;   % volatility of variance
params.rho = rho;   % correlation between Brownian motions (NOTE: methods which assume rho=0 will display in output)
model = 6;
modelInput = getModelInput(model, T/M, r, q, params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analytical (Bernard and Cui, 2014)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/STOCHASTIC_VOL/DV_Swaps_Options/Analytical_Swaps')
tic
[ref, KcH] = hestonfairstrike(r, v0, theta, eta, Sigmav, T, rho, M);
time_ana = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Time-Changed Markov Chain Approximation (Cui, Kirkby, Nguyen, 2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/STOCHASTIC_VOL/Helper_Functions')
addpath('../../PROJ/TIME_CHANGED/DV_Swaps_Options')

% ==========================
% CTMC Params
% ==========================
ParamsCtmc.varGridMult = .05;
ParamsCtmc.gamma = 6;  % Heston gamma = 4 is good for T ~ 1
ParamsCtmc.Nx = 70; %the number of Markov states
n = 0; % number of time steps in time disretization... set to 0 to do continuous time version

ParamsDiffus.model = 1;
ParamsDiffus.eta = eta; ParamsDiffus.rho = rho; ParamsDiffus.theta = theta; ParamsDiffus.Sigmav = Sigmav; ParamsDiffus.v0 = v0;

hFunc = @(u) u;   % tau = int h(X_s) ds
levyExponent = @(z) -0.5*1i*z - 0.5*z.^2; 
tic;
price_tc_mc = TC_Levy_VarianceSwap(r, q, T, M, levyExponent, hFunc, n, ParamsDiffus, ParamsCtmc);
time_tc_mc = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SV-PROJ  (Cui, Kirkby, Nguyen, 2017)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/STOCHASTIC_VOL/DV_Swaps_Options/')
addpath('../../PROJ/STOCHASTIC_VOL/Helper_Functions/')
N  = 2^9;    % Number of basis elements
L1 = 14;     % Truncation width parameter
m_0 = 30;      % Number of variance states in CTMC
gamma = 5.5;     % State space grid width parameter  (5.5 is a decent starting point, but depends on the model)
varGridMult = .8;    % State space stretching paramter, chose in (0,1)... closer to 1 is more uniform grid
gridMethod = 4;     % Determines the grid method... 4 is a good choice


jumpParams = {}; jumpParams.Nothing = 0; psi_J = @(u)0*[u>0];

alph = GetAlph_DisreteVariance( 0, 0, 1, params, T, L1 );
tic;
price_sv_ctmc = PROJ_DiscreteVariance_StochVol( N,alph,M,r,T,0,m_0,psi_J, 1, params, gridMethod, gamma, varGridMult, 1 );
time_sv_ctmc = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Monte Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Monte_Carlo/')
addpath('../../Monte_Carlo/Variance_Swaps_Options/')

tic
N_sim = 10^4; disc = exp(-r*T); scheme = 5;  S_0 = 100; % NOTE: S_0 is irrelevant
mult = 2;
Spath = Simulate_Heston_Euler_Schemes( N_sim, M*mult, T, S_0, r, q, params, scheme);
[price_MC, stdErr] = Price_MC_Var_Swaps_func(Spath, disc, M, mult);
price_MC_L = price_MC - 2*stdErr; price_MC_U = price_MC + 2*stdErr;
time_MC = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n---------------------------------------------\n')
fprintf('Method       | Price       |    Err   |  CPU  \n')
fprintf('---------------------------------------------\n')
fprintf('Analytical   | %.8f  | %.2e | %.4f  \n', ref, 0, time_ana)
fprintf('TC-MC (rho=0)| %.8f  | %.2e | %.4f  \n', price_tc_mc, abs(price_tc_mc-ref), time_tc_mc)
fprintf('SV-PROJ      | %.8f  | %.2e | %.4f  \n', price_sv_ctmc, abs(price_sv_ctmc-ref), time_sv_ctmc)
fprintf('MC-Euler     |[%.3f,%.3f]| %.2e | %.4f \n', price_MC_L, price_MC_U, abs(price_MC-ref), time_MC)


