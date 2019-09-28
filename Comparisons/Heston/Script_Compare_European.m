%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EUROPEAN OPTION PRICE COMPARISON (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Compare Methods For European Options Under Hestons Model
% Author:      Justin Kirkby
% 
% Methods: 1) Kahl-Jackel-Lord Approach
%          2) PROJ
%          3) Monte Carlo, Using Lord et al (2010) Schemes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../../PROJ/LEVY/European_Options')
addpath('../../PROJ/LEVY/RN_CHF')
addpath('../../PROJ/LEVY/Helper_Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call = 1;    %For call use 1 (2 for put)
S_0  = 100;  %Initial price
W    = 100;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = .00;  %Interest rate (NOTE: set to zero for comparison with Kahl-Jackel-Lord, based on Forward price)
q    = .00;  %dividend yield (NOTE: keep this at zero for now)
T    = 1;    %Time (in years)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS  (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = {};

params.v0 = 0.0175; % initial variance
params.theta = 0.0398;   % long term variance level
params.eta = 1.5768;   % rate of variance mean reversion
params.Sigmav = 0.5751;   % volatility of variance
params.rho = -0.5711;   % correlation between Brownian motions

modelInput = getModelInput(6, T, r, q, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Kahl-Jackel-Lord Approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Fourier/Heston/')
tic
price_KJL = Heston1993KahlJaeckelLordRev3(call, S_0,W,T,0,r,q, params.v0, params.theta, params.rho, params.eta, params.Sigmav);
time_KJL = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROJ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order = 3;  %Choose spline order from { 0,1,2,3} => {Haar, Linear, Quadratic, Cubic}
logN  = 12;   %Uses N = 2^logN  gridpoint 
L1 = 18;

% ----------------------
N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]
alpha = getTruncationAlpha(T, L1, modelInput, 6);

tic
price_PROJ = PROJ_European(order, N, alpha, r, q, T, S_0, W, call, modelInput.rnCHF, modelInput.c1*T);
time_PROJ = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Monte Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../Monte_Carlo/')
addpath('../../Monte_Carlo/European/')

tic
N_sim = 10^5; M = 800; disc = exp(-r*T); scheme = 5;
Spath = Simulate_Heston_Euler_Schemes( N_sim, M, T, S_0, r, q, params, scheme);
[price_MC, stdErr] = Price_MC_European_Strikes_func(Spath, disc, call, W );
price_MC_L = price_MC - 2*stdErr; price_MC_U = price_MC + 2*stdErr;
time_MC = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----------------------------\n')
fprintf('Method   | Price       | CPU  \n')
fprintf('-----------------------------\n')
fprintf('PROJ     | %.8f  | %.4f  \n', price_PROJ, time_PROJ)
fprintf('KJL      | %.8f  | %.4f  \n', price_KJL, time_KJL)
fprintf('MC-Euler |[%.3f,%.3f]| %.4f \n', price_MC_L, price_MC_U, time_MC)


