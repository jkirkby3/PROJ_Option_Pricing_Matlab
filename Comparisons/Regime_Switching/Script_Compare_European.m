%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Script to Compare Methods For European Options Under Regime Switching Model
%              Compares accuracy/CPU of the following methods,
%                   1) Monte Carlo
%                   2) PROJ
%                      ... More to come
%
% Author:      Justin Kirkby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);

% ---------------------
%  Contract/Market Params
% ---------------------
call = 1;    %For call use 1 (else, its a put)
S_0  = 100;  %Initial price
r    = .05;  %Interest rate
q    = .00;  %dividend yield
T    = 1;    %Time (in years)
W    = 100;  %strike to price

% ---------------------
% Regime Switching Diffusion Params
% ---------------------
% Transition Matrix (dictates how the regimes transition)
Q = [-1 0.5 0.5;
    0.5 -1 0.5; 
    0.5 0.5 -1];  

drift_vec = [r-q  r-q  r-q];  % Drift in each state
sigma_vec = [0.15  0.25  0.35]; % Volatility in each state

initial_state = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Monte Carlo Pricer (Euler)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Variance reduction techniques should be used in practice
addpath('../../Monte_Carlo/')
addpath('../../Monte_Carlo/European')
N_sim = 10^5;
M = 800;

Spath = Simulate_RegimeSwitching_Diffusion_func( N_sim, M, T, S_0, drift_vec, sigma_vec, Q, initial_state);

disc = exp(-r*T);
[price_MC, stdErr] = Price_MC_European_Strikes_func(Spath, disc, call, W );
price_MC_L = price_MC - 2*stdErr; price_MC_U = price_MC + 2*stdErr;
time_MC = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  PROJ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/REGIME_SWITCHING/')
addpath('../../PROJ/REGIME_SWITCHING/European_Options')
logN  = 8;   %Uses N = 2^logN  gridpoint 
L1 = 8;  % determines grid witdth (usually set L1 = 8 to 15 for Levy)

alpha = L1*sqrt(T)*max(sigma_vec);   % Choose grid width based on the largest volatility
N = 2^logN;    % grid roughly centered on [- alph, alph]

tic
price_PROJ = PROJ_RegimeSwitching_European(3, N, alpha, r, q, T, S_0, W, call, Q, sigma_vec, initial_state);
time_PROJ = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n---------------------------------------------\n')
fprintf('Method      | Price           |  CPU \n')
fprintf('---------------------------------------------\n')
fprintf('PROJ        | %.8f     | %.4f \n', price_PROJ,  time_PROJ)
fprintf('MC-Euler    | [%.3f,%.3f] | %.4f \n', price_MC_L, price_MC_U, time_MC)

