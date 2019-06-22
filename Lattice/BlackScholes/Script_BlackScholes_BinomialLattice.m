%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BINOMIAL LATTICE OPTION PRICER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price European/American options in Black-Scholes Models
%              using Binomial Lattice
% Author:      Justin Kirkby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K     = 100;    % strike
S_0   = 100;    % initial stock price
r     = 0.05;   % interest rate
T     = 1;      % time to maturity
sigma = 0.15;   % volatility of diffusion
M     = 252;    % Number of time steps
call  = 0;      % call (set to 1) else put
american = 1;   % American (set to 1) else European

price_binomial = BinomialLattice_BlackScholes_func( S_0,K,r,T,sigma,M, call, american)
price_trinomial = TrinomialLattice_BlackScholes_func( S_0,K,r,T,sigma,M, call, american)

    
    