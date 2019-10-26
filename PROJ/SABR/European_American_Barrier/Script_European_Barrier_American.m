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
addpath('../Helper_Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -----------------
% Contract/Market Params
% -----------------
T    = 1;
r    = 0.00;
S_0  = 1.1;   % Futures price 
Kvec = S_0*[0.6 0.8 0.90 0.95 1.00 1.05 1.10 1.2 1.4];
call = 0;    % 1 = call option, else put
M    = 100;   % number of time steps

contract_type = 1;  % 1 = European, 2 = American, 3 = Down and Out Barrier
L    = 0.6*S_0;     % For barrier contract, this is the barrier

% -----------------
% Numerical Params
% -----------------
CTMCParams.m_0 = 30;
CTMCParams.N = 90;
CTMCParams.gridMult_v = 0.5;
CTMCParams.gridMult_s = 0.05;  %Grid mult param for S 
CTMCParams.gamma = 6;  %Grid width param for variance grid

% -----------------
% Model Params
% -----------------
ModParams.beta   = .6;
ModParams.alpha  = 0.08;
ModParams.v0     = 0.2;
ModParams.rho    = 0;

% ModParams.beta   = .6;
% ModParams.alpha  = 0.3;
% ModParams.v0     = 0.25;
% ModParams.rho    = -0.5;

% ModParams.beta   = .7;
% ModParams.alpha  = 0.08;
% ModParams.v0     = 0.2;
% ModParams.rho    = 0;

% ModParams.beta   = .3;
% ModParams.alpha  = 0.6;
% ModParams.v0     = 0.4;
% ModParams.rho    = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
prices = SABR_EurBarAmer_func(call, M, T, S_0, Kvec, r, CTMCParams, ModParams, contract_type, L);
toc
fprintf('%.8f \n', prices)

plot(Kvec / S_0, prices)
ylabel('price', 'interpreter', 'latex')
xlabel('moneyness, $K/S_0$', 'interpreter', 'latex')
grid on;

