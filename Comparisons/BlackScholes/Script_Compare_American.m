%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AMERICAN OPTION PRICE COMPARISON (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Compare Methods For American (Bermudan) Options Under Black Scholes Model
%              This script compares accuracy/CPU of the following methods,
%                   Lattices (Binomial / Trinomial)
%                   Monte Carlo (Longstaff-Schwartz)
%                   PROJ
%
% Author:      Justin Kirkby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../../PROJ/LEVY/RN_CHF')
addpath('../../PROJ/LEVY/Helper_Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS (Price PUT Option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_0  = 100;  %Initial price
W    = 100;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = .05;  %Interest rate
T    = 1;    %Time (in years)
M    = 1000;  %Num exercise dates for bermudan

sigma = 0.15;  % volatility of diffusion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1; params = {}; params.sigmaBSM = sigma; q = 0;
modelInput = getModelInput(model, T/M, r, q, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  PROJ Method (use PROJ for high accuracy reference)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/LEVY/American_Options')

call = 0; % NOTE: this script only supports PUT option
logN  = 10;   %Uses N = 2^logN  gridpoint 
L1 = 12;
alpha = getTruncationAlpha(T, L1, modelInput, model);

% ----------------------
ref = PROJ_Bermudan_Put(M, S_0, W, r, T, modelInput.rnCHF, 2^16, alpha);
% ----------------------

tic
N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]
price_PROJ = PROJ_Bermudan_Put(M, S_0, W, r, T, modelInput.rnCHF, N, alpha);
time_PROJ = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Binomial / Trinomial Lattice Pricer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Lattice/BlackScholes/')
M_steps = M;  % number of binomial time steps
tic
price_binom = BinomialLattice_BlackScholes_func(S_0, W, r, T, sigma, M_steps, call, 1);
time_binom = toc;

tic
price_trinom = TrinomialLattice_BlackScholes_func(S_0, W, r, T, sigma, M_steps, call, 1);
time_trinom = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Monte Carlo - Lonstaff Schwartz (MC-LS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Monte_Carlo/')
addpath('../../Monte_Carlo/American')
N_sim = 10^4;
Spath = Simulate_Jump_Diffusion_func( N_sim, M, T, S_0, r, q, sigma);

dt = T/M; disc = exp(-r*dt);
[price_MC, stdErr] = Price_MC_American_Strikes_func(Spath, disc, call, W );

price_MC_L = price_MC - 2*stdErr; price_MC_U = price_MC + 2*stdErr;
time_MC = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n---------------------------------------------\n')
fprintf('Method      |    Price    |    Err   |  CPU \n')
fprintf('---------------------------------------------\n')
fprintf('Ref         | %.8f  |          |       \n', ref)
fprintf('---------------------------------------------\n')
fprintf('PROJ        | %.8f  | %.2e | %.4f \n', price_PROJ, abs(ref-price_PROJ), time_PROJ)
fprintf('Binomial    | %.8f  | %.2e | %.4f \n', price_binom, abs(ref-price_binom), time_binom)
fprintf('Trinomial   | %.8f  | %.2e | %.4f \n', price_trinom, abs(ref-price_trinom), time_trinom)
fprintf('---------------------------------------------\n')
fprintf('MC-LS       |[%.3f,%.3f]| %.2e | %.4f \n', price_MC_L, price_MC_U, abs(ref-price_MC), time_MC)
fprintf('---------------------------------------------\n')
