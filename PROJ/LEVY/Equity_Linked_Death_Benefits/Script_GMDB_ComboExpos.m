%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GMWB PRICER (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Guaranteed Minimum Death Benefits (GMDB)
%              using the PROJ method, under combination of exponential mortaily model
%
% Author:      Zhimin Zhang  (Original Code)
%              Justin Lars Kirkby (Convert into common framework)
%
% References:  (1) Valuing Equity-Linked Death Benefits in General Exponential
%               Levy Models, J. Comput. and Appl. Math. 2019 (Z. Zhang, Y. Yong, W. Yu)
%              (2) Efficient Option Pricing By Frame Duality with The Fast
%               Fourier Transform, SIAM J. Financial Math., 2015 (J.L. Kirkby)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../RN_CHF')
addpath('../Helper_Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call = 0;    %For call use 1 (else, its a put)
S_0  = 100;  %Initial price
W    = 120;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = .05;  %Interest rate
T    = 20;    %Time (in years), set to -1 to price a perpetual (ie no time to expiry)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS  (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1;   %See Models Below (e.g. model 1 is Black Scholes), and choose specific params

params = {};
params.model = model;

if model == 1 %BSM (Black Scholes Merton)
    params.sigmaBSM = 0.25;    %CHOOSE   

elseif model == 3 %NIG
    params.sigma  = 0.25;
    params.alpha = 2;
    params.beta  = 0.5;
    params.delta = 0.05;
    
elseif model == 4 %MJD (Merton Jump Diffusion)
    params.sigma  = 0.25;
    params.lam    = 0.6;
    params.muj    = 0.01;
    params.sigmaj = 0.13;
    
elseif model == 5 %Kou Double Expo
    params.sigma = 0.25;
    params.lam   = 0.6;
    params.p_up  = 0.5;
    params.eta1  = 1;
    params.eta2  = 4;
    
elseif model == 8 % Variance Gamma
    params.sigmaGBM = 0.25; % geometric brownian motion add-on
    params.theta = 0.01;
    params.sigma = 0.05;
    params.nu = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 3) CHOOSE Time-Until-Death PARAMETERS (Combination of expos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params_death = {};
params_death.A = [3 -2];  % coefficients of combination of expos
params_death.lambda = [0.08 0.12];  % expo rates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 4) CHOOSE PROJ PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order = 1;  %Spline order (for now just linear)
P = 8;
Pbar = P+2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if order == 1
    tic
    price = PROJ_GMDB_ComboExpos_Linear(P, Pbar, S_0, W, call, r, params, params_death, T);
    toc
    fprintf('%.8f \n', price)
else
    fprintf('Only orders 1 (linear) and 2 (quadratic) are supported')
end




