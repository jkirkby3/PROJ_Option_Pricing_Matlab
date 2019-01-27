[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);


addpath('../RN_CHF')
addpath('../Helper_Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cliquet/Equity Index Annuity (EIA) PRICER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Barrier options in Levy Models
%              using the PROJ method
% Author:      Justin Kirkby
% References:  (1) Equity-linked Annuity pricing with Cliquet-style
%               guarantees in regime-switching and stochastic volatility
%               models with jumps, Insurance: Math. and Economics, 2017
%              (2) Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform, SIAM J. Financial Math., 2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r    = .05;  %Interest rate
q    = .00;  %dividend yield
T    = 1;    %Time (in years)
M    = 12;  %number of discrete monitoring points

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contract: 1 = sum of local caps
%           2 = sum of local caps & floors
%           3 = cliquet: local & global caps & floors
%           4 = cliquet: local floor & cap, global floor, NO GLOBAL CAP  
%           5 = MPP: ie monthly point-to-point or Monthly Cap Sum (Bernard, Li)
%           6 = Multiplicative Cliquet (e.g. Hieber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

contract = 3;  % 3 is the standard contract type (the general Cliquet)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
contractParams.K   = 1;  % "Strike" (principal)

if contract ~= 6
    contractParams.C  = .06;  % Local cap
    contractParams.CG = 0.75*M*contractParams.C;  % Global Cap
    contractParams.F  = 0.01;   % Local Floor
    contractParams.FG = 1.25*M*contractParams.F;  % Global Floor
else
    % Mutliplicative style cliquet
    contractParams.Alpha = .25; % exponent of return
    contractParams.C  = 1.05;
    contractParams.F  = 1.00; %NOTE: must be greater than zero
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 1;   %See Models Below (e.g. model 1 is Black Scholes), and choose specific params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CHOOSE PROJ PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UseCumulant = 1;  %Set to 1 to use the cumulant base rule (Approach 1) to determine gridwidth, else used fixed witdth (Approach 2)

%---------------------
% APPROACH 1: Cumulant Based approach for grid width
% (see "Robust Option Pricing with Characteritics Functions and the BSpline Order of Density Projection")
%---------------------
if UseCumulant ==1  %With cumulant based rule, choose N and Alpha (N = 2^(P+Pbar) based on second approach)
    logN  = 11;   %Uses N = 2^logN  gridpoint 
    L1 = 12;  % determines grid witdth (usually set L1 = 8 to 15 for Levy, or 18 for Heston)
%---------------------
% APPROACH 2: Manual GridWidth approach 
%--------------------- 
else %Manually specify resolution and Pbar
    P     = 7;  % resolution is 2^P
    Pbar  = 3;  % Determines density truncation grid with, 2^Pbar 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CHOOSE MODEL PARAMETERS 
%%%  Note: rnCHF is the risk netural CHF, c1,c2,c4 are the cumulants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelInput = getModelInput(model, T/M, r, q, params);

if UseCumulant ==1  % Choose density truncation width based on cumulants
    alpha = getTruncationAlpha(T, L1, modelInput, model);
else    % Manually supply density truncation width above
    logN = P + Pbar;
    alpha = 2^Pbar/2;
end
N = 2^logN;    

tic
if contract == 6
    price = MultiplicativeCliquet_PROJ( N,alpha,M,r,T,modelInput.rnCHF, contract,contractParams);
else
    price = Cliquet_LEVY_PROJ(N,alpha,M,r,q,T,modelInput.rnCHF,contract,contractParams);

end
toc

fprintf('price: %.8f\n', price)


