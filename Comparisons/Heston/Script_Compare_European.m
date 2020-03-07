%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EUROPEAN OPTION PRICE COMPARISON (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Compare Methods For European Options Under Hestons Model
% Author:      Justin Kirkby
% 
% Methods: 1) Kahl-Jackel-Lord Approach (Fourier, Heston Model)
%          2) PROJ (Kirkby, 2015), European Levy/Heston 
%          3) CONV (Lord, Fang, Bervoets, Oosterlee, 2008), European Levy/Heston 
%          4) Carr-Madan (2008), European Levy/Heston Pricer 
%          5) Regime Switching Fourier PROJ (Cui, Kirkby, Nguyen, 2017) - Stoch Vol / RS Pricer
%          6) Time-Changed Markov Chain (Cui, Kirkby, Nguyen, 2019), assumes rho = 0
%          7) Monte Carlo, Using Lord et al (2010) Low-Bias Schemes
%               ... More to come
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../../PROJ/LEVY/European_Options')
addpath('../../PROJ/LEVY/RN_CHF')
addpath('../../PROJ/LEVY/Helper_Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call = 1;    %For call use 1 (2 for put)
S_0  = 100;  %Initial price
W    = 100;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = .00;  %Interest rate (NOTE: set to zero for comparison with Kahl-Jackel-Lord, based on Forward price)
q    = .00;  %dividend yield (NOTE: keep this at zero for now)
T    = 0.5;    %Time (in years)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS  (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = {};

params.v0 = 0.02; % initial variance
params.theta = 0.02;   % long term variance level
params.eta = 1.6;   % rate of variance mean reversion
params.Sigmav = 0.3;   % volatility of variance
params.rho = 0;   % correlation between Brownian motions (NOTE: methods which assume rho=0 will display in output)

modelInput = getModelInput(6, T, r, q, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Kahl-Jackel-Lord (KJL) Approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Fourier/Heston/')
tic
price_KJL = Heston1993KahlJaeckelLordRev3(call, S_0,W,T,0,r,q, params.v0, params.theta, params.rho, params.eta, params.Sigmav);
time_KJL = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROJ (Kirkby, 2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logN  = 12;   %Uses N = 2^logN  gridpoint 
L1 = 30;

% ----------------------
N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]
alpha = getTruncationAlpha(T, L1, modelInput, 6);

tic
price_PROJ = PROJ_European(3, N, alpha, r, q, T, S_0, W, call, modelInput.rnCHF, modelInput.c1*T);
time_PROJ = toc;

ref = PROJ_European(3, 2^15, 2*alpha, r, q, T, S_0, W, call, modelInput.rnCHF, modelInput.c1*T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Time-Changed Markov Chain Approximation (Cui, Kirkby, Nguyen, 2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/TIME_CHANGED/European/')
addpath('../../PROJ/STOCHASTIC_VOL/Helper_Functions')
params.model = 1;

ParamsCtmc.varGridMult = .01;
ParamsCtmc.gamma = 6;  % Heston gamma = 3 is good for T ~ 1
ParamsCtmc.Nx = 100; %the number of Markov states

ProjParams.order = 3;
ProjParams.alph = 2^2;
ProjParams.N_proj = 2^8;

n = 0; % number of time steps in time disretization... set to 0 to do continuous time version
hFunc = @(u) u;   % tau = int h(X_s) ds
levyExponent = @(z) -0.5*1i*z - 0.5*z.^2; 

tic
price_TCMC = PROJ_TimeChanged_Levy_European(r,q,S_0,T,W,call, levyExponent, hFunc, n, params, ParamsCtmc, ProjParams);
time_TCMC = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SV-PROJ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/STOCHASTIC_VOL/European/')
addpath('../../PROJ/STOCHASTIC_VOL/Barrier/')
addpath('../../PROJ/STOCHASTIC_VOL/Helper_Functions/')
% This version uses the Stoch Vol pricer for Barrier options to price European (More of a multiple purpose method)
if call == 1
    down = 1; H = S_0 / 8;   % TODO: this needs to account for the variance of the underlying.
else
    down = 0; H = S_0 * 8;
end

N             = 2^10;    %number of points in density expansion... Value grid size is K:=N/2
m_0           = 40;  % number of CTMC grid points
gamma         = 5;  % CTMC grid width param
gridMethod    = 4; gridMultParam = 0.2; M = 1; psi_J = @(u)0*[u>0];
alpha         = 5;

tic
price_SVP = Barrier_StochasticVol_func(N,alpha,call,down,S_0,W,H,M,r,T,m_0,psi_J,1, params, gridMethod, gamma, gridMultParam);
time_SVP = toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Carr-Madan Fourier Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Fourier/CarrMadan/')
N = 2^15;
tic
price_CM = CarrMadan_European_Price_Strikes(S_0, W, modelInput.rnCHF, N, T, r, q, call);
time_CM = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CONV Fourier Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Fourier/CONV/')
N = 2^16;
tic
price_CONV = CONV_European_Price(S_0, W, modelInput.rnCHF, T, r, call, N, alpha);
time_CONV = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Monte Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Monte_Carlo/')
addpath('../../Monte_Carlo/European/')

tic
N_sim = 10^5; M = 800; disc = exp(-r*T); scheme = 5;
Spath = Simulate_Heston_Euler_Schemes( N_sim, M, T, S_0, r, q, params, scheme);
[price_MC, stdErr] = Price_MC_European_Strikes_func(Spath, disc, call, W );
price_MC_L = price_MC - 2*stdErr; price_MC_U = price_MC + 2*stdErr;
time_MC = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n---------------------------------------------\n')
fprintf('Method       | Price       |    Err   |  CPU  \n')
fprintf('---------------------------------------------\n')
fprintf('PROJ         | %.8f  | %.2e | %.4f  \n', price_PROJ, abs(price_PROJ-ref), time_PROJ)
fprintf('KJL          | %.8f  | %.2e | %.4f  \n', price_KJL, abs(price_KJL-ref), time_KJL)
fprintf('CONV         | %.8f  | %.2e | %.4f  \n', price_CONV, abs(price_CONV-ref), time_CONV)
fprintf('Carr-Madan   | %.8f  | %.2e | %.4f  \n', price_CM, abs(price_CM-ref), time_CM)
fprintf('TC-MC (rho=0)| %.8f  | %.2e | %.4f  \n', price_TCMC, abs(price_TCMC-ref), time_TCMC)
fprintf('SV-PROJ      | %.8f  | %.2e | %.4f  \n', price_SVP, abs(price_SVP-ref), time_SVP)
fprintf('MC-Euler     |[%.3f,%.3f]| %.2e | %.4f \n', price_MC_L, price_MC_U, abs(price_MC-ref), time_MC)


