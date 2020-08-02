%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2D EUROPEAN/Barrier/Bermudan OPTION PRICER (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price European options in 2D Diffusion Models
%              using the CTMC Approximation
%
%              This script is configured initially to compare against analytical exchange option price,
%              but it can be used to price numerous payoffs with European/Barrier/Bermudan style
%
% Author:      Justin Kirkby
% References:  (1) A General Continuous Time Markov Chain Approximation for
%               Multi-Asset option pricing with systems of correlated diffusions,
%               Applied Math. and Comput., 2020 (with Duy Nguyen and Dang Nguyen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../')
addpath('../../Analytical/BlackScholes/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0.0;               % interest rate
S_0s = [97 97];        % initial asset prices
sigmas = [0.15 0.15];  % volatilities per asset
qs = [0 0];            % div yeilds per asset
rho = 0.5;             % correlation between brownian motions 
T = 1;                 % time to maturity
M = 10;                % num monitoring points


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 3) Choose Contract Params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contractParams = {};

% ---------------
% Contract Types:
% ---------------  
% 1: European, Single Step Pricing
% 2: European, Multi Step Pricing (M above controls number of steps)
% 3: Bermudan (M above controls number of monitoring points)
% 4: Barrier (M above controls number of monitoring points)
contractParams.contract = 1;

if contractParams.contract == 4
    % Set Barrier For Barrier Option (initially configured to price with only a barrier on the first underlying)
    contractParams.barriers_1 = [0 50];
    contractParams.barriers_2 = [0 50000000000];
end

% ---------------
% Payoff Types:
% ---------------  
% 1: Linear, G = S_1  (linear payoff in first underlying)
% 2: Linear, G = S_2  (linear payoff in second underlying)
% 3: Exchange, G = (S_1 - S_2)^+
% 4: Spread,  G = (S_1 - S_2 - K)^+   (NOTE: must set strike, K)
% 5: Geometric Basket Call / Put,  G = (sqrt(S_1) * sqrt(S_2) - K)^+  (for the call)
% 6: Arithmetic Basket Call / Put,  G = (sqrt(S_1) * sqrt(S_2) - K)^+  (for the call)
% 7: Call-on-Max and Put-on-Min, Gcall = (max(S_1,S_2) - K)^+ , Gput = (K - min(S_1,S_2))^+
% 8: Call/put on just S_2, G = (S_2 - K)^+  (for the call)
% 9: Best-of / Worst-of,  G = max(S_1,S_2), G = min(S_1,S_2)
contractParams.payoff_type = 4; 

% Payoff Dependent Contract parameters (ignored for certain payoffs)
contractParams.K = 0;    % Strike(depending on payoff type)       
contractParams.call = 1;  % call/put

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 4) CHOOSE CTMC PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = {};
params.num_devs = 5;  % num std devs used in the grid
params.GridMultParam = 0.1;  % non-uniformity param, in (0,1)
params.gridMethod = 7;  % choose the grid method (several RnD versions)
params.m_0 = 120;  % num CTCM states



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 5) PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[vals, c_index_1, c_index_2, y_1, y_2] = price_2d_ctmc( S_0s, T, r, rho, sigmas, qs, params, contractParams, M);
price_ctmc = vals(c_index_1, c_index_2);
time_ctmc = toc;

fprintf('CTMC Price: %.8f   (time: %.3f)\n', price_ctmc, time_ctmc)

% Analytical price in case of Spread 
if contractParams.payoff_type == 4 && (contractParams.contract <= 2)  &&  contractParams.K == 0
    price_exchange = Price_Exchange_Option_Margrabe_2D( S_0s(1), S_0s(2), T, rho, sigmas(1), sigmas(2), qs(1), qs(2));
    fprintf('Exact Exchange Price: %.8f\n', price_exchange)
    fprintf('Error: %.2e\n', price_ctmc - price_exchange)
end

