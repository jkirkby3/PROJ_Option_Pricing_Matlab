%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2D Exchange OPTION PRICE COMPARISON (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to compare pricing methods for Exchange options under 2D Black-Schoels
%              This script compares accuracy/CPU of the following methods,
%                   CTMC approximation
%                   Monte Carlo
%                   Magrabe (anlytical/exact)
%                   COS Method
%
%
% Author:      Justin Kirkby
% References:  (1) A General Continuous Time Markov Chain Approximation for
%                 Multi-Asset option pricing with systems of correlated diffusions,
%                 Applied Math. and Comput., 2020 (JL Kirkby, Duy Nguyen, Dang Nguyen)
%              (2) M. J. Ruijter and C. W. Oosterlee. Two-dimensional Fourier cosine
%                 series expansion method for pricing financial options. SIAM Journal on
%                 Scientific Computing, 34(5):B642–B671, 2012.
%              (3) Margrabe, W. (1978): "The value of an option to exchange one asset for another,"
%                 Journal of Finance, 33 (1978), pp. 177-186.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../../Analytical/BlackScholes/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0.05;              % interest rate
S_0s = [97 97];        % initial asset prices
sigmas = [0.15 0.15];  % volatilities per asset
rho = 0.5;             % correlation between brownian motions 
T = 1;                 % time to maturity
call = 1;              % call = 1, else put


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Magrabe Analytical in case of Exchange option 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qs = [0 0];            % div yeilds per asset (NOTE: not all methods support this, keep at 0)

ref = Price_Exchange_Option_Margrabe_2D( S_0s(1), S_0s(2), T, rho, sigmas(1), sigmas(2), qs(1), qs(2));


fprintf('\n--------------------------------------------------\n')
fprintf('Method       |    Price    |    Err   |  CPU \n')
fprintf('--------------------------------------------------\n')
fprintf('Ref (Magrabe)| %.8f  |          |       \n', ref)
fprintf('--------------------------------------------------\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COS (code available for K==0 only, ie exchange option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../../Fourier/COS')
tic
N = 22; L = 9;
price_cos = Price_Exchange_COS2D(S_0s,T,r,sigmas(1),sigmas(2),rho,N,L);
time_cos = toc;

fprintf('COS-2D       | %.8f  | %.2e | %.4f \n', price_cos, abs(ref-price_cos), time_cos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CTMC Approximation Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../CTMC/')
addpath('../../CTMC/Diffusion_2D')
contractParams = {};
contractParams.contract = 1;      % 1: European
contractParams.payoff_type = 4;   % 4: Spread,  G = (S_1 - S_2 - K)^+   (NOTE: must set strike, K)
contractParams.K = 0;            
contractParams.call = call;       
M = 10;                         % num monitoring points

params = {};
params.num_devs = 5;  % num std devs used in the grid
params.GridMultParam = 0.1;  % non-uniformity param, in (0,1)
params.gridMethod = 7;  % choose the grid method (several RnD versions)
params.m_0 = 120;  % num CTCM states

tic
[vals, c_index_1, c_index_2, y_1, y_2] = price_2d_ctmc( S_0s, T, r, rho, sigmas, qs, params, contractParams, M);
price_ctmc = vals(c_index_1, c_index_2);
time_ctmc = toc;

fprintf('CTMC         | %.8f  | %.2e | %.4f \n', price_ctmc, abs(ref-price_ctmc), time_ctmc)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Monte Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Monte_Carlo')

tic
N_sim = 10^5; M = 200; 
dt = T / M;
drifts = [r r] - qs;
exponential = 1;  % see code, this means we are simulating geometric brownian motion
[paths_1, paths_2] = Simulate_Diffusion_2D(S_0s, drifts, sigmas, rho,  N_sim, M, dt, exponential);
S_1 = paths_1(:,M + 1);
S_2 = paths_2(:,M + 1);

payoffs = max(S_1 - S_2, 0);
price_MC = exp(-r*T)*mean(payoffs);
stderr = exp(-r*T)*std(payoffs) / sqrt(N_sim);
price_MC_L = price_MC - 2*stderr; price_MC_U = price_MC + 2*stderr;
time_MC = toc;

fprintf('MC           |[%.3f,%.3f]| %.2e | %.4f \n', price_MC_L, price_MC_U, abs(ref-price_MC), time_MC)
fprintf('--------------------------------------------------\n')

