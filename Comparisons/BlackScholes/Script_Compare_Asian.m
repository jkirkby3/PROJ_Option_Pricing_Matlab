%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ASIAN OPTION PRICE COMPARISON (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Compare Methods For Arithmetc Asian Options Under Black Scholes Model
%              This script compares accuracy/CPU of the following methods,
%                   PROJ
%                   Monte Carlo
%                   Vorst Approx
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
M    = 100;  %Num exercise dates for bermudan
call = 1;    %call = 1, put = -1

sigma = 0.15;  % volatility of diffusion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1; params = {}; params.sigmaBSM = sigma; q = 0;
modelInput = getModelInput(model, T/M, r, q, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Vorst Approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Approx/BlackScholes/Asian')
tic
price_vorst = Price_Asian_VorstApprox_BlackScholes(S_0, sigma, M, W, call, T, r, q);
time_vorst = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  PROJ Method (use PROJ for high accuracy reference)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/LEVY/Asian_Options')

logN = 7;   %Uses N = 2^logN  gridpoint 
L1 = 8;
alpha = getTruncationAlpha(T, L1, modelInput, model);

% ----------------------
ref = PROJ_Asian(2^10, alpha, S_0, M, W, call, T, r, q, modelInput.rnCHF, modelInput.RNmu*T/M);
% ----------------------

tic
N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]
price_PROJ = PROJ_Asian(N, alpha, S_0, M, W, call, T, r, q, modelInput.rnCHF, modelInput.RNmu*T/M);
time_PROJ = toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Monte Carlo (Euler)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Monte_Carlo/')
addpath('../../Monte_Carlo/Asian')
N_sim = 10^6;
tic
Spath = Simulate_Jump_Diffusion_func( N_sim, M, T, S_0, r, q, sigma);

dt = T/M; disc = exp(-r*T);
[price_MC, stdErr] = Price_MC_Asian_Strikes_func(Spath, call, disc, W, M, 1);

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
fprintf('Vorst       | %.8f  | %.2e | %.4f \n', price_vorst, abs(ref-price_vorst), time_vorst)
fprintf('---------------------------------------------\n')
fprintf('MC          |[%.3f,%.3f]| %.2e | %.4f \n', price_MC_L, price_MC_U, abs(ref-price_MC), time_MC)
fprintf('---------------------------------------------\n')
