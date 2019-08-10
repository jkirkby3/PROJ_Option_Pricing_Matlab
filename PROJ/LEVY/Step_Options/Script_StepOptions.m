%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP OPTION PRICER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Step (Soft) Barrier options in Levy Models
%              using the PROJ method
%  Payoff is exp(-stepRho*tau_M)*(S_T - W)^+ for a call, where tau_M is the amount of time spent in knock-out region
% Author:      Justin Kirkby
% References:   
%              (1) Robust Barrier Option Pricing by Frame Projection under
%              Exponential Levy Dynamics, App. Math. Finance, 2017
%              (2) Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform, SIAM J. Financial Math., 2015
%              (3) Step Options, (Linetsky, V.), Math. Finance 1999.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../RN_CHF')
addpath('../Helper_Functions')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1): CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_0  = 100;  %Initial price
W    = 100;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = .05;  %Interest rate
q    = .0;   %Dividend yield
T    = 1;    %Time (in years)
call = 0;    %For call use 1 (else, its a put)
down = 1;    %down-out or up-out (down=1 => down-and-out)
H    = 90;   %barrier
M    = 52;  %number of discrete monitoring points 

stepRho = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2): CHOOSE MODEL PARAMETERS (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1;   %See Models Below (e.g. model 1 is Black Scholes), and choose specific params
params = {};

if model == 1 %BSM (Black Scholes Merton)
    params.sigmaBSM = 0.18;    %CHOOSE   
    
elseif model == 2 %CGMY
    params.C  = 0.02; 
    params.G  = 5; 
    params.MM = 15; 
    params.Y  = 1.2;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 3) CHOOSE PROJ PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logN  = 10;   %Uses N = 2^logN  gridpoint 
L1 = 10;  % determines grid witdth (usually set L1 = 8 to 15 for Levy, or 18 for Heston)

% For Automated Parameter adjustment
alphMult = 1.1;
TOLProb = 5e-08;
TOLMean = 1e-05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Note: rnCHF is the risk netural CHF, c1,c2,c4 are the cumulants
modelInput = getModelInput(model, T/M, r, 0, params);
N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]

c1 = modelInput.c1; c2 = modelInput.c2; c4 = modelInput.c4; rnCHF = modelInput.rnCHF; rnCHF_T = modelInput.rnCHF_T;


if (down == 1 && call ~= 1) || (down ~=1 && call == 1)
    tic
    price = PROJ_StepOption_AutoParam(N,stepRho,call,down, S_0,W,H,M,r,q,rnCHF,T,L1,c2,c4, alphMult,TOLProb,TOLMean,rnCHF_T);
    fprintf('PRICE: %.8f \n', price)
    toc
else
    fprintf('Sorry, currently only Up and out calls, and down and out puts have been coded \n')
end


