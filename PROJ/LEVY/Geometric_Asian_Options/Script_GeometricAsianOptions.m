%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GEOMETRIC ASIAN OPTION PRICER  (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Geometric Asian options in Levy Models using the PROJ method
% Author:      Justin Kirkby
%
% Reference:   (1) An Efficient Transform Method For Asian Option Pricing, SIAM J. Financial Math., 2016
%              (2) Efficient Option Pricing by Frame Duality with the Fast Fourier Transform. 
%                  SIAM J. Financial Math (2015), Kirkby, J.L
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../RN_CHF')
addpath('../Helper_Functions')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call = 1;    %For call use 1 (else, its a put)
S_0  = 100;  %Initial price
W    = 100;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = 0.05;  %Interest rate
q    = 0.00;  %dividend yield
T    = 1;    %Time (in years)
M    = 52;   %Number of monitoring points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1;   %See Models Below (e.g. model 1 is Black Scholes), and choose specific params
params = {};

if model == 1 %BSM (Black Scholes Merton)
    params.sigmaBSM = 0.15;    %CHOOSE   
    
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
    
elseif model == 8 % Variance Gamma 
    params.sigma = 0.2; 
    params.nu = 0.85;  
    params.theta = 0;    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 3) CHOOSE PROJ PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UseCumulant = 1;  %Set to 1 to use the cumulant base rule (Approach 1) to determine gridwidth, else used fixed witdth (Approach 2)

%---------------------
% APPROACH 1: Cumulant Based approach for grid width
% (see "Robust Option Pricing with Characteritics Functions and the BSpline Order of Density Projection")
%---------------------
if UseCumulant ==1  %With cumulant based rule, choose N and Alpha (N = 2^(P+Pbar) based on second approach)
    logN  = 8;   %Uses N = 2^logN  gridpoint 
    L1 = 10;  % determines grid witdth (usually set L1 = 8 to 15 for Levy)
%---------------------
% APPROACH 2: Manual GridWidth approach 
%--------------------- 
else %Manually specify resolution and Pbar
    P     = 6;  % resolution is 2^P
    Pbar  = 3;  % Determines density truncation grid with, 2^Pbar 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Note: rnCHF is the risk netural CHF, c1,c2,c4 are the cumulants
dt = T/M;
modelInput = getModelInput(model, dt, r, q, params);

if UseCumulant ==1  % Choose density truncation width based on cumulants
    alpha = getTruncationAlpha(T, L1, modelInput, model);
else    % Manually supply density truncation width above
    logN = P + Pbar;
    alpha = 2^Pbar/2;
end
N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]

tic
price = PROJ_Geometric_Asian(N, alpha, S_0, M, W, call, T, r, q, modelInput.rnSYMB);
toc


fprintf('Geometric Price: %.8f \n', price)

compare_arithmetic = 1;

if compare_arithmetic
    addpath('../Asian_Options')
    price_arith = PROJ_Asian( N,alpha,S_0,M,W,call,T,r,q, modelInput.rnCHF, modelInput.RNmu*dt);
    fprintf('Arithmetic Price: %.8f \n', price_arith)
end
    
