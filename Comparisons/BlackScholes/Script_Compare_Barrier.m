%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% KNOCK-OUT BARRIER OPTION PRICE COMPARISON (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Compare Methods For Barrier Options Under Black Scholes Model
%              This script compares accuracy/CPU of the following methods,
%                   Monte Carlo
%                   PROJ
%
% Author:      Justin Kirkby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../../PROJ/LEVY/RN_CHF')
addpath('../../PROJ/LEVY/Helper_Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_0    = 100;   %Initial price
W      = 100;   %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r      = 0.05;  %Interest rate
q      = 0.02;  %dividend yield
T      = 1;     %Time (in years)
call   = 1;     %For call use 1 (else, its a put)
down   = 1;     %down-out or up-out (down=1 => down-and-out)
H      = 90;    %barrier (Knock-Out)
M      = 52;    %number of discrete monitoring points
rebate = 0;     % rebate paid immediately upon passing the barrier (knocking-out) 

sigma = 0.15;  % volatility of diffusion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1; params = {}; params.sigmaBSM = sigma;
modelInput = getModelInput(model, T/M, r, q, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Reference Price (Using PROJ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/LEVY/Barrier_Options')
L1_ref = 20; N_ref = 2^17;
alpha_ref = getTruncationAlpha(T, L1_ref, modelInput, model);
price_ref = PROJ_Barrier(N_ref, alpha_ref, call, down, S_0, W, H, M, r, q, modelInput.rnCHF, T, rebate); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  PROJ Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 2^12;    % grid roughly centered on [c1 - alph, c1 + alph]
L1 = 10;
alpha = getTruncationAlpha(T, L1, modelInput, model);

tic
price_PROJ = PROJ_Barrier(N, alpha, call, down, S_0, W, H, M, r, q, modelInput.rnCHF, T, rebate); 
time_PROJ = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Monte Carlo Pricer (Exact Sim, Multi-Step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Variance reduction techniques should be used in practice
addpath('../../Monte_Carlo/')
addpath('../../Monte_Carlo/Barrier')
N_sim = 5*10^5;  % number of paths
mult = 3;   % multiplier for simulation (see below) to reduce bias

tic
M_mult = M*mult;  %time partitioning to reduce bias
Spath = Simulate_Jump_Diffusion_func( N_sim, M_mult + 1, T, S_0, r, q, sigma);

[price_MC, stdErr] = Price_MC_Barrier_Strikes_func(Spath, call, down, H, W, M, mult, rebate, r, T);

price_MC_L = price_MC - 2*stdErr; price_MC_U = price_MC + 2*stdErr;
time_MC = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n---------------------------------------------\n')
fprintf('Method      |    Price    |    Err   |  CPU \n')
fprintf('---------------------------------------------\n')
fprintf('Reference   | %.8f  |          |       \n', price_ref)
fprintf('---------------------------------------------\n')
fprintf('PROJ        | %.8f  | %.2e | %.4f \n', price_PROJ, abs(price_ref-price_PROJ), time_PROJ)
fprintf('---------------------------------------------\n')
fprintf('MC-Exact    |[%.3f,%.3f]| %.2e | %.4f \n', price_MC_L, price_MC_U, abs(price_ref-price_MC), time_MC)
fprintf('---------------------------------------------\n')

