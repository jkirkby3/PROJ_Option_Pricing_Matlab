%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EUROPEAN OPTION PRICE COMPARISON (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Compare Methods For European Options Under Black Scholes Model
%              This script compares accuracy/CPU of the following methods,
%                   Monte Carlo (Euler, Exact Sim)
%                   Lattices (Binomial, Trinomial)
%                   PROJ
%                   Finite Difference (Explicit, Implicit, Crank-Nicholson)
%                   Fourier (Carr-Madan, CONV)
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
%%%  Analytical Pricer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Analytical/BlackScholes/')
price_True = BSM_Greeks(0, S_0, sigma, r, q, T, W, call);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Binomial / Trinomial Lattice Pricer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Lattice/BlackScholes/')
M = 1500;  % number of binomial time steps
tic
price_binom = BinomialLattice_BlackScholes_func(S_0, W, r, T, sigma, M, call, 0);
time_binom = toc;

tic
price_trinom = TrinomialLattice_BlackScholes_func(S_0, W, r, T, sigma, M, call, 0);
time_trinom = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Finite Difference (Explicit / Implicit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PDE_FiniteDifference/BlackScholes/')
Smax = 5*S_0; Smin = 0;

% Explicit
tic
dS = 1;  dt = dS^2 / Smax^2 / sigma^2;  % Ensure stability
price_explicitFD = ExplicitFD_BlackScholes_func(S_0, W, r, T, sigma, call, dS, dt, Smax, Smin);
time_explicitFD = toc;

% Fully Implicit
tic
dS = 0.5; dt = 1/200;
price_implicitFD = ImplicitFD_BlackScholes_func(S_0, W, r, T, sigma, call, dS, dt, Smax, Smin);
time_implicitFD = toc;

% Crank Nicolson (Implicit/Explicit)
tic
dS = 0.5; dt = 1/200;
price_crankNicFD = CrankNicolsonFD_BlackScholes_func(S_0, W, r, T, sigma, call, dS, dt, Smax, Smin);
time_crankNicFD = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Carr-Madan Fourier Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Fourier/CarrMadan/')
N = 2^15;
tic
price_CM = CarrMadan_European_Price_Strikes(S_0, W, modelInput.rnCHF, N, T, r, q, call);
time_CM = toc;

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
price_PROJ = PROJ_European(order, N, alpha, r, q, T, S_0, W, call, modelInput.rnCHF, modelInput.c1*T);
time_PROJ = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CONV Fourier Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Fourier/CONV/')
N = 2^16;
tic
price_CONV = CONV_European_Price(S_0, W, modelInput.rnCHF, T, r, call, N, alpha);
time_CONV = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Monte Carlo Pricer (Euler)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Variance reduction techniques should be used in practice
addpath('../../Monte_Carlo/')
addpath('../../Monte_Carlo/European')
N_sim = 10^5;  % number of paths
M = 500; % number of time steps for euler
disc = exp(-r*T);  % dicount factor
tic
Spath = Simulate_Jump_Diffusion_func( N_sim, M, T, S_0, r, q, sigma);
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
fprintf('Method      |    Price    |    Err   |  CPU \n')
fprintf('---------------------------------------------\n')
fprintf('Exact       | %.8f  |          |       \n', price_True)
fprintf('---------------------------------------------\n')
fprintf('PROJ        | %.8f  | %.2e | %.4f \n', price_PROJ, abs(price_True-price_PROJ), time_PROJ)
fprintf('CONV        | %.8f  | %.2e | %.4f \n', price_CONV, abs(price_True-price_CONV), time_CONV)
fprintf('Carr-Madan  | %.8f  | %.2e | %.4f \n', price_CM, abs(price_True-price_CM), time_CM)
fprintf('---------------------------------------------\n')
fprintf('Binomial    | %.8f  | %.2e | %.4f \n', price_binom, abs(price_True-price_binom), time_binom)
fprintf('Trinomial   | %.8f  | %.2e | %.4f \n', price_trinom, abs(price_True-price_trinom), time_trinom)
fprintf('---------------------------------------------\n')
fprintf('FD-Explicit | %.8f  | %.2e | %.4f \n', price_explicitFD, abs(price_True-price_explicitFD), time_explicitFD)
fprintf('FD-Implicit | %.8f  | %.2e | %.4f \n', price_implicitFD, abs(price_True-price_implicitFD), time_implicitFD)
fprintf('FD-CrankNic | %.8f  | %.2e | %.4f \n', price_crankNicFD, abs(price_True-price_crankNicFD), time_crankNicFD)
fprintf('---------------------------------------------\n')
fprintf('MC-Euler    |[%.3f,%.3f]| %.2e | %.4f \n', price_MC_L, price_MC_U, abs(price_True-price_MC), time_MC)
fprintf('MC-Exact    |[%.3f,%.3f]| %.2e | %.4f \n', price_EMC_L, price_EMC_U, abs(price_True-price_EMC), time_EMC)
fprintf('---------------------------------------------\n')

