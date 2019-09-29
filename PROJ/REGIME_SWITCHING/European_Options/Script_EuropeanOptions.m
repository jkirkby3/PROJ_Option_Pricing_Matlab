%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EUROPEAN OPTION PRICER (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price European options in Regime Switching Diffusion Models
%              using the PROJ method
% Author:      Justin Kirkby
% References:  (1) Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform, SIAM J. Financial Math., 2015
%              (2) A unified approach to Bermudan and Barrier options under stochastic
%               volatility models with jumps. J. Economic Dynamics and Control, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call = 1;    %For call use 1 (else, its a put)
S_0  = 100;  %Initial price
W    = 100;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = .05;  %Interest rate
q    = .00;  %dividend yield
T    = 1;    %Time (in years)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transition Matrix (dictates how the regimes transition)
Q = [-1 0.5 0.5;
    0.5 -1 0.5; 
    0.5 0.5 -1];  

sigma_vec = [0.15  0.25  0.35]; % Volatility in each state

initial_state = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 3) CHOOSE PROJ PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order = 3;  %Choose spline order from { 0,1,2,3} => {Haar, Linear, Quadratic, Cubic}, Haar is least accurate
logN  = 8;   %Uses N = 2^logN  gridpoint 
L1 = 8;  % determines grid witdth (usually set L1 = 8 to 15 for Levy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = L1*sqrt(T)*max(sigma_vec);   % Choose grid width based on the largest volatility
N = 2^logN;    % grid roughly centered on [- alph, alph]

tic
price = PROJ_RegimeSwitching_European(order, N, alpha, r, q, T, S_0, W, call, Q, sigma_vec, initial_state);
toc

fprintf('%.8f \n', price)

