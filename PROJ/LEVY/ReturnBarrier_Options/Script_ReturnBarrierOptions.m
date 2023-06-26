%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Return BARRIER OPTION PRICER  (This script demomstrates convergence of
%%%                                the pricing algorithm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Script to Price Return Barrier options in Levy Models
%              using the PROJ method
% Author:      Justin Kirkby
% Reference: (1) "The Return Barrier and Return Timer Option with Pricing Under
%              Levy Processes", J.L. Kirkby and J-P Aguilar, Expert Systems with Applications, 2023
%            (2) "Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform", J.L. Kirkby, SIAM J. Financial Math., 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../RN_CHF')
addpath('../Helper_Functions')
addpath('../European_Options')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r    = 0.03; % Interest rate
q    = 0.00; % dividend yield
T    = 1;    % Time (in years)
S_0  = 100;  % Initial underlying (spot) value
W    = 100;  % Strike
call = 0;    % For call use 1 (else, its a put)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daily Monitoring Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = -0.09;   % lower return barrier
u = 0.09;    % upper return barrier
M = 252;     % number of discrete monitoring points

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monthly Monitoring Example
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% l = -0.2;  % lower return barrier
% u = 0.2;   % upper return barrier
% M = 12;    % number of discrete monitoring points

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Weekly Monitoring Example
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% l = -0.15;  % lower return barrier
% u = 0.15;   % upper return barrier
% M = 52;     % number of discrete monitoring points

contractParams.F   = l;
contractParams.C   = u;

contractParams.K   = W;  % strike
contractParams.call = call;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1;   %See Models Below (e.g. model 1 is Black Scholes), and choose specific params

% BSM, MJD, NIG, GGMY, KDE, VG
%  1    4    3     2    5    8
case_id = 1;
params = getModel(model,case_id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 3) CHOOSE PROJ PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logN  = 11;   % Uses N = 2^logN  gridpoint 
L1_T = 16;    % determines grid witdth (usually set L1 = 8 to 15 for Levy)
L1_dt = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Note: rnCHF is the risk netural CHF, c1,c2,c4 are the cumulants
dt = T/M;
modelInput = getModelInput(model, dt, r, q, params, T);

alpha_T = getTruncationAlpha(T, L1_T, modelInput, model);
alpha_dt = getTruncationAlpha(dt, L1_dt, modelInput, model);

N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]
order = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% European Price for comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
price_eur = PROJ_European(order, N, alpha_T, r, q, T, S_0, W, call, modelInput.rnCHF_T, modelInput.c1);
fprintf('European Price: %.8f \n', price_eur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
price_ref = PROJ_ReturnBarrier(2^13, alpha_T, alpha_dt, M, r, q, T, S_0,  modelInput.rnCHF, contractParams);
fprintf('Ref RBO Price: %.8f \n', price_ref)

Nvec = 2.^[6 7 8 9 10 11];

vals_cubic = zeros(1, length(Nvec));
errs_cubic = zeros(1, length(Nvec));
time_cubic = zeros(1, length(Nvec));

N_trials = 1;

for i=1:length(vals_cubic)
   tic;
   for j = 1:N_trials
     val = PROJ_ReturnBarrier(Nvec(i), alpha_T, alpha_dt, M, r, q, T, S_0,  modelInput.rnCHF, contractParams);
   end
   time = toc;
   time_cubic(i) = time/ N_trials;
   
   vals_cubic(i) = val;
   errs_cubic(i) = abs(val-price_ref);
end

rates_cubic = -log2(abs(errs_cubic(2:end)./errs_cubic(1:end-1)));

for i=1:length(vals_cubic)
    if i > 1
        fprintf('%.0f &  %.10f & %.2e & % .2f & %.5f \\\\ \n', log2(Nvec(i)), vals_cubic(i), errs_cubic(i), rates_cubic(i-1), time_cubic(i));
    else
        fprintf('%.0f &  %.10f & %.2e & -- & %.5f \\\\ \n', log2(Nvec(i)), vals_cubic(i), errs_cubic(i), time_cubic(i));
    end
end

