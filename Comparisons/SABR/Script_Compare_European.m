%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SABR European option pricer comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Compare European Option pricers in SABR Model
% Author:      Justin Kirkby
% References:  (1) General Valuation Framework for SABR and Stochastic Local Volatility
%                   Models. SIAM J. Financial Mathematics, 2018.
%              (2) Full Fledge SABR with Markov Chains. Wimott, 2019.
% Disclaimer: this is research code, not production code. The parameter settings etc 
%            should not be expected to work as is in all cases. Determining the right parameter
%            settings will depend on your usecase. This is left up to the user. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T    = 1;         % time to matruity
r    = 0.05;      % interest rate
F_0  = 1.1;       % Futures price 
K    = F_0*0.95;  % strike
call = 1;         % 1 = call option, else put

% -----------------
% Model Params
% -----------------
ModParams.beta   = .6;
ModParams.alpha  = 0.08;
ModParams.v0     = 0.2;
ModParams.rho    = 0.0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Reference Price (Using Obloj Asymptotic Formula with correction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Approx/SABR')
addpath('../../Analytical/BlackScholes')

tic
IV = SABR_Hagan_Obloj_ImpliedVol( K, F_0, ModParams.v0, T, ModParams.alpha, ModParams.beta, ModParams.rho );
price_obloj = BSM_Greeks( 0, F_0/exp(r*T), IV, r, 0, T, K, call)
time_oblog=toc;

price_ref = price_obloj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Antonov Approximation (Using Asymptotic Formula)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
price_ant = SABR_European_AntonovApprox(F_0,K,T,call,r,ModParams)
time_ant = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Monte Carlo (Euler)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Monte_Carlo')
addpath('../../Monte_Carlo/European')
M = 400;
N_sim = 5*10^4;

tic;
Spath = Simulate_SLV_func( N_sim, M, T, F_0, 0, 0, 2, ModParams);
disc = exp(-r*T);
[price_MC, stdErr] = Price_MC_European_Strikes_func(Spath, disc, call, K );

price_MC_L = price_MC - 2*stdErr; price_MC_U = price_MC + 2*stdErr;
time_MC = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CTMC-Double Layer Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/SABR/European_American_Barrier/')
addpath('../../PROJ/SABR/Helper_Functions/')

% -----------------
% Numerical Params
% -----------------
M    = 100;   % number of time steps
contract_type = 1;    % 1 = European, 2 = American, 3 = Down and Out Barrier
L    = 0.6*F_0;     % For barrier contract, this is the barrier

CTMCParams.m_0 = 26;
CTMCParams.N = 90;
CTMCParams.gridMult_v = 0.5;
CTMCParams.gridMult_s = 0.05;  %Grid mult param for S 
CTMCParams.gamma = 6;  %Grid width param for variance grid

tic
price_ctmc = SABR_EurBarAmer_func(call, M, T, F_0, K, r, CTMCParams, ModParams, contract_type, L)
time_ctmc = toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n---------------------------------------------\n')
fprintf('Method      |    Price    |    Err   |  CPU \n')
fprintf('---------------------------------------------\n')
fprintf('Ref (Obloj) | %.8f  |          |       \n', price_ref)
fprintf('---------------------------------------------\n')
fprintf('Antonov     | %.8f  | %.2e | %.4f \n', price_ant, abs(price_ref-price_ant), time_ant)
fprintf('CTMC        | %.8f  | %.2e | %.4f \n', price_ctmc, abs(price_ref-price_ctmc), time_ctmc)
fprintf('MC-Euler    |[%.3f,%.3f]| %.2e | %.4f \n', price_MC_L, price_MC_U, abs(price_ref-price_MC), time_MC)
fprintf('---------------------------------------------\n')




