%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% American/Bermudan Option Pricier for SABR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Bermudan Options in SABR Model using CTMC Method
% Author:      Justin Kirkby
% References:  (1) General Valuation Framework for SABR and Stochastic Local Volatility
%                   Models. SIAM J. Financial Mathematics, 2018.
% Disclaimer: this is research code, not production code. The parameter settings etc 
%            should not be expected to work as is in all cases. Determining the right parameter
%            settings will depend on your usecase. This is left up to the user. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
%addpath('../../STOCHASTIC_VOL/Helper_Functions')
addpath('../Helper_Functions')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -----------------
% Contract/Market Params
% -----------------
T    = 1;
r    = 0.00;
S0   = 1.1;
Kvec = [1.0 1.06 1.1 1.12];
L = 0.6;     % For barrier contract, this is the barrier
call = 0;    % 1 = call option, else put
M = 100;   % number of time steps
contract_type = 1;  % 1 = European, 2 = American, 3 = Down and Out Barrier


% -----------------
% Numerical Params
% -----------------
CTMCParams.m_0 = 30;
CTMCParams.N = 90;
CTMCParams.gridMult_v = 0.5;
CTMCParams.gridMult_s = 0.05;  %Grid mult param for S 
CTMCParams.gamma = 4;  %Grid width param for variance grid

% -----------------
% Model Params
% -----------------
ModParams.beta   = .7;
ModParams.alpha  = 0.08;
ModParams.v0     = 0.2;
ModParams.rho    = -0.4;

% ModParams.beta   = .3;
% ModParams.alpha  = 0.6;
% ModParams.v0     = 0.4;
% ModParams.rho    = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
prices = SABR_EurBarAmer_func(call, M, T, S0, Kvec, r, CTMCParams, ModParams, contract_type, L);
toc
fprintf('%.8f \n', prices)


