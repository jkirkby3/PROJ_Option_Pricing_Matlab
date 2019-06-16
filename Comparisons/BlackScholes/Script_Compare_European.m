%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EUROPEAN OPTION PRICE COMPARISON (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Compare Methods For European Options Under Black Scholes Model
% Author:      Justin Kirkby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../../PROJ/LEVY/RN_CHF')
addpath('../../PROJ/LEVY/Helper_Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call = 1;    %For call use 1 (else, its a put)
S_0  = 100;  %Initial price
W    = 100;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = .05;  %Interest rate
T    = 1;    %Time (in years)

sigma = 0.15;  % volatility of diffusion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1; params = {}; params.sigmaBSM = sigma; q = 0;
modelInput = getModelInput(model, T, r, q, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  PROJ Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/LEVY/European_Options')
order = 3;  %Choose spline order from { 0,1,2,3} => {Haar, Linear, Quadratic, Cubic}
logN  = 11;   %Uses N = 2^logN  gridpoint 
L1 = 16;

% ----------------------
N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]
alpha = getTruncationAlpha(T, L1, modelInput, model);

tic
price_PROJ = PROJ_European( order,N,alpha,r,q,T,S_0,W,call, modelInput.rnCHF, modelInput.c1*T);
time_PROJ = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Analytical Pricer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Analytical/BlackScholes/')
price_True = BSM_Greeks( 0, S_0, sigma, r, q, T, W, call);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Lattice Pricer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Lattice/BlackScholes/')
M = 2000;  % number of binomial time steps
tic
price_binom = BinomialLattice_BlackScholes_func( S_0, W, r, T, sigma, M, call, 0);
time_binom = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Monte Carlo Pricer (Euler)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Variance reduction techniques should be used in practice
addpath('../../Monte_Carlo/')
addpath('../../Monte_Carlo/European')
jumpModel = 0; jumpParams = {};
N_sim = 10^5;  % number of paths
M = 500; % number of time steps for euler
disc = exp(-r*T);  % dicount factor
tic
Spath = Simulate_Jump_Diffusion_func( N_sim, M, T, S_0, r, q, sigma, jumpModel, jumpParams);
[price_MC, stdErr] = Price_MC_European_Strikes_func(Spath, disc, call, W );
price_MC_L = price_MC - 2*stdErr; price_MC_U = price_MC + 2*stdErr;
time_MC = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Monte Carlo - EXACT Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the trivial case of European option under Black Scholes, we can EXACTLY simulated the terminal distriution
% NOTE: Variance reduction techniques should be used in practice
N_sim = 5*10^5;  % number of paths
tic
W1 = randn(N_sim,1); 
Spath = S_0*exp((r - q - .5*sigma^2)*T +  sigma*sqrt(T)*W1);   % Exact simulation
[price_EMC, stdErr] = Price_MC_European_Strikes_func(Spath, disc, call, W );
price_EMC_L = price_EMC - 2*stdErr; price_EMC_U = price_EMC + 2*stdErr;
time_EMC = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n---------------------------------------------\n')
fprintf('Method   |    Price    |    Err   |  CPU \n')
fprintf('---------------------------------------------\n')
fprintf('Exact    | %.8f  |          |       \n', price_True)
fprintf('PROJ     | %.8f  | %.2e | %.4f \n', price_PROJ, abs(price_True-price_PROJ), time_PROJ)
fprintf('Lattice  | %.8f  | %.2e | %.4f \n', price_binom, abs(price_True-price_binom), time_binom)
fprintf('MC-Euler |[%.3f,%.3f]| %.2e | %.4f \n', price_MC_L, price_MC_U, abs(price_True-price_MC), time_MC)
fprintf('MC-Exact |[%.3f,%.3f]| %.2e | %.4f \n', price_EMC_L, price_EMC_U, abs(price_True-price_EMC), time_EMC)


