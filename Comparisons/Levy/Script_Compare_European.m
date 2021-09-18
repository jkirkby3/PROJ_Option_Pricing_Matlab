%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EUROPEAN OPTION PRICE COMPARISON (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Compare Methods For European Options Under Levy Models
%              This script compares accuracy/CPU of the following Fourier methods:
%                   1) PROJ (2015)
%                   2) Carr-Madan (2000)
%                   3) CONV (2008)
%                   4) Lewis (2001) / Lipton (2002)
%                   5) Mellin Transform (Aguilar, 2019)
%                       ... More to come
%   
% Author:      Justin Kirkby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../../PROJ/LEVY/RN_CHF')
addpath('../../PROJ/LEVY/Helper_Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call = 0;    %For call use 1 (else, its a put)
S_0  = 4000;  %Initial price
W    = 4000;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = .01;  %Interest rate
q    = .00;  %Dividend yield
T    = 1;    %Time (in years)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS  (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 8;   %See Models Below (e.g. model 1 is Black Scholes), and choose specific params

params = {};
if model == 1 %BSM (Black Scholes Merton)
    params.sigmaBSM = 0.15;    %CHOOSE   
    
elseif model == 2 %CGMY
    params.C  = 0.02; 
    params.G  = 5; 
    params.MM = 15; 
    params.Y  = 1.2;

elseif model == 3 %NIG
    params.alpha = 8.9932;
    params.beta  = -4.5176;
    params.delta = 1.1528;
    
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
    
elseif model == 6 % Heston Model  
    params.v_0 = 0.0175; % initial variance
    params.theta = 0.0398;   % long term variance level
    params.kappa =1.5768;   % rate of variance mean reversion
    params.sigma_v = 0.5751;   % volatility of variance
    params.rho = -0.5711;   % correlation between Brownian motions
    
elseif model == 7 %KoBoL
    params.c  = 0.02; 
    params.lam_p  = 15; 
    params.lam_m = -5; 
    params.nu  = 1.2;
    
elseif model == 8 % Variance Gamma 
    params.sigma = 0.2; 
    params.nu = 0.85;  
    params.theta = -0.1;   

elseif model == 9 % Bilateral Gamma 
    params.alpha_p = 1.18; 
    params.lam_p = 10.57;  
    params.alpha_m = 1.44; 
    params.lam_m = 5.57;
end

modelInput = getModelInput(model, T, r, q, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Carr-Madan Fourier Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Fourier/CarrMadan/')
N = 2^16;
tic
price_CM = CarrMadan_European_Price_Strikes(S_0, W, modelInput.rnCHF, N, T, r, q, call);
time_CM = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  PROJ Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../PROJ/LEVY/European_Options')
logN  = 12;   %Uses N = 2^logN  gridpoint 
if model == 6
    L1 = 28;
elseif model == 1
    L1 = 18;
else
    L1 = 16;
end


% ----------------------
N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]
alpha = getTruncationAlpha(T, L1, modelInput, model);

tic
price_PROJ = PROJ_European(3, N, alpha, r, q, T, S_0, W, call, modelInput.rnCHF, modelInput.c1*T);
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
%%%  Hilbert Transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Fourier/HilbertTransform/')
N = 2^14; Nh = 2^5;
h = 2*pi/Nh;

tic
price_Hilb = Hilbert_European_Price(h, N, r, q, T, S_0, W, call, modelInput.rnCHF);
time_Hilb = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Lewis (2001) Fourier Method  (see also Lipton (2002))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Fourier/Lewis/')
tic
price_Lewis = Lewis_European_Price(S_0, W, modelInput.rnCHF, T, r, q, call, 500);
time_Lewis = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Mellin Transform (J-P Aguilar, 2020)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Fourier/MellinTransform/')
has_mellin = 0;
if model == 3  % NIG
    tic
    N_terms = 50;
    price_Mellin = Mellin_NIG_European_Price( S_0, W, T, r, q, call, params.alpha, params.beta, params.delta, N_terms);
    time_Mellin = toc;
    has_mellin = 1;
end
if model == 8 %  Variance Gamma
    tic
    N_terms = 24;
    price_Mellin = Mellin_VG_European_Price( S_0, W, T, r, q, call, params.sigma, params.theta, params.nu, N_terms);
    time_Mellin = toc;
    has_mellin = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Reference Price (Using PROJ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nref = 2^16; L1Ref = max(16, L1+4);   % Params to obtain reference price
alphaRef = getTruncationAlpha(T, L1Ref, modelInput, model);
price_Ref = PROJ_European(3, Nref, alphaRef, r, q, T, S_0, W, call, modelInput.rnCHF, modelInput.c1*T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n---------------------------------------------\n')
fprintf('Method      |    Price      |    Err   |  CPU \n')
fprintf('---------------------------------------------\n')
fprintf('Reference   | %.8f  |          |       \n', price_Ref)
fprintf('---------------------------------------------\n')
fprintf('PROJ        | %.8f  | %.2e | %.4f \n', price_PROJ, abs(price_Ref-price_PROJ), time_PROJ)
fprintf('CONV        | %.8f  | %.2e | %.4f \n', price_CONV, abs(price_Ref-price_CONV), time_CONV)
fprintf('Carr-Madan  | %.8f  | %.2e | %.4f \n', price_CM, abs(price_Ref-price_CM), time_CM)
fprintf('Lewis(2001) | %.8f  | %.2e | %.4f \n', price_Lewis, abs(price_Ref-price_Lewis), time_Lewis)
fprintf('Hilbert     | %.8f  | %.2e | %.4f \n', price_Hilb, abs(price_Ref-price_Hilb), time_Hilb)
if has_mellin == 1
   fprintf('Mellin      | %.8f  | %.2e | %.4f \n', price_Mellin, abs(price_Ref-price_Mellin), time_Mellin) 
end
fprintf('---------------------------------------------\n')

