%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOOKBACK / HINDSIGHT OPTION PRICER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Lookback/Hindsight options in Levy Models
%              using the PROJ method
% Author:      Justin Kirkby
% References:  (1) American and Exotic Option Pricing with Jump Diffusions
%               and other Levy Processes, J. Computational Finance 2018
%              (2) Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform, SIAM J. Financial Math., 2015
%
[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../RN_CHF')
addpath('../Helper_Functions')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Floating Strike (Lookback) Put:  max{S_m: 0<=m<=M} - S_T
% Floating Strike (Lookback) Call: S_T - min{S_m: 0<=m<=M}
% Fixed Strike (Hindsight) Put: (W - min{S_m: 0<=m<=M})^+
% Fixed Strike (Hindsight) Call: (max{S_m: 0<=m<=M} - W)^+

S_0  = 100;  %Initial price
W    = 100;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = .1;   %Interest rate
q    = .05;  %dividend yield
T    = 0.5;  %Time (in years)

M    = 50;

%%%%%%%%%%%%%%%%                      
% NOTE: currently Only floating strike Put, and fixed Strike Call Are Supported
% The other two contract types will be finished soon (hopefully)
%%%%%%%%%%%%%%%%
call = 1;
floating_strike = 0;  % used 1 for floating stike, otherwise is fixed strike
                      % for floating strike, W (strike) is irrelevant param
if ~(floating_strike == 1 && call ~= 1) && ~(floating_strike ~= 1 && call == 1)
    error("This contract type not currently finished. Implementaion in progress");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1;   %See Models Below (e.g. model 1 is Black Scholes), and choose specific params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CHOOSE PROJ PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order = 3;  %Choose spline order from { 0,1,2,3} => {Haar, Linear, Quadratic, Cubic}
UseCumulant = 1;  %Set to 1 to use the cumulant base rule (Approach 1) to determine gridwidth, else used fixed witdth (Approach 2)

%---------------------
% APPROACH 1: Cumulant Based approach for grid width
% (see "Robust Option Pricing with Characteritics Functions and the BSpline Order of Density Projection")
%---------------------
if UseCumulant ==1  %With cumulant based rule, choose N and Alpha (N = 2^(P+Pbar) based on second approach)
    logN  = 13;   %Uses N = 2^logN  gridpoint 
    L1 = 12;  % determines grid witdth (usually set L1 = 8 to 15 for Levy, or 18 for Heston)
%---------------------
% APPROACH 2: Manual GridWidth approach 
%--------------------- 
else %Manually specify resolution and Pbar
    P     = 7;  % resolution is 2^P
    Pbar  = 2;  %Determines grid with, 2^Pbar 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CHOOSE MODEL PARAMETERS 
%%%  Note: rnCHF is the risk netural CHF, c1,c2,c4 are the cumulants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = {};

if model == 1 %BSM (Black Scholes Merton)
    params.sigmaBSM = 0.3;    %CHOOSE   
    
elseif model == 2 %CGMY
    params.C  = 4; 
    params.G  = 50; 
    params.MM = 60; 
    params.Y  = 0.7;

elseif model == 3 %NIG
    params.alpha = 15;
    params.beta  = -5;
    params.delta = 0.5;
    
elseif model == 4 %MJD (Merton Jump Diffusion)
    params.sigma  = 0.12;
    params.lam    = 0.4;
    params.muj    = -0.12;
    params.sigmaj = 0.18;
    
elseif model == 5 %Kou Double Expo
    params.sigma = 0.15;
    params.lam   = 3;
    params.p_up  = 0.2;
    params.eta1  = 25;
    params.eta2  = 10;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelInput = getModelInput(model, T/M, r, q, params);

if UseCumulant ==1  % Choose density truncation width based on cumulants
    alpha = getTruncationAlpha(T, L1, modelInput, model);
else    % Manually supply density truncation width
    logN = P + Pbar;
    alpha = 2^Pbar/2;
end

N = 2^logN;

tic
price = PROJ_Lookback( N,alpha, S_0,W,call,r,q,M,T,modelInput.rnSYMB, floating_strike); 
toc


fprintf('%.8f \n', price)
