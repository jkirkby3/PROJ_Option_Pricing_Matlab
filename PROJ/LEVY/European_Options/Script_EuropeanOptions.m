%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EUROPEAN OPTION PRICER (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price European options in Levy/Heston Models
%              using the PROJ method
% Author:      Justin Kirkby
% References:  (1) Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform, SIAM J. Financial Math., 2015
%              (2) Robust Option Pricing with Characteristic Functions and
%              the B-Spline Order of density Projection, JCF, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../RN_CHF')
addpath('../Helper_Functions')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call = 1;    %For call use 1 (else, its a put)
S_0  = 100;  %Initial price
W    = 100;  %Strike            %NOTE: no error handling in place for extreme values of W (increase grid if strike falls outside)
r    = .05;  %Interest rate
q    = .01;  %dividend yield
T    = 1;    %Time (in years)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS  (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
elseif model == 6 % Heston Model  
    params.v_0 = 0.0175; % initial variance
    params.theta = 0.0398;   % long term variance level
    params.kappa =1.5768;   % rate of variance mean reversion
    params.sigma_v = 0.5751;   % volatility of variance
    params.rho = -0.5711;   % correlation between Brownian motions
    
elseif model == 8 % Variance Gamma 
    params.sigma = 0.2; 
    params.nu = 0.85;  
    params.theta = 0.1; 
    
elseif model == 9 % Bilateral Gamma 
    params.alpha_p = 1.18; 
    params.lam_p = 10.57;  
    params.alpha_m = 1.44; 
    params.lam_m = 5.57;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 3) CHOOSE PROJ PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order = 3;  %Choose spline order from { 0,1,2,3} => {Haar, Linear, Quadratic, Cubic}
UseCumulant = 1;  %Set to 1 to use the cumulant base rule (Approach 1) to determine gridwidth, else used fixed witdth (Approach 2)

%---------------------
% APPROACH 1: Cumulant Based approach for grid width
% (see "Robust Option Pricing with Characteritics Functions and the BSpline Order of Density Projection")
%---------------------
if UseCumulant ==1  %With cumulant based rule, choose N and Alpha (N = 2^(P+Pbar) based on second approach)
    logN  = 14;   %Uses N = 2^logN  gridpoint 
    if model == 6  % Heston
        L1 = 18;
    else
        L1 = 12;  % determines grid witdth (usually set L1 = 8 to 15 for Levy, or 18 for Heston)
    end
%---------------------
% APPROACH 2: Manual GridWidth approach 
%--------------------- 
else %Manually specify resolution and Pbar
    P     = 7;  % resolution is 2^P
    Pbar  = 3;  % Determines density truncation grid with, 2^Pbar 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Note: rnCHF is the risk netural CHF, c1,c2,c4 are the cumulants
modelInput = getModelInput(model, T, r, q, params);

if UseCumulant ==1  % Choose density truncation width based on cumulants
    alpha = getTruncationAlpha(T, L1, modelInput, model);
else    % Manually supply density truncation width above
    logN = P + Pbar;
    alpha = 2^Pbar/2;
end
N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]


tic
price = PROJ_European( order,N,alpha,r,q,T,S_0,W ,call, modelInput.rnCHF, modelInput.c1*T);
toc

fprintf('%.8f \n', price)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_convergence = 1;
plot_smile = 1;

if plot_convergence
    figure();
    Nvec = 2.^[3 4 5 6 7 8 9 10];
    errs = zeros(1, length(Nvec));
   
    ref = PROJ_European( order,2^12,alpha,r,q,T,S_0,W ,call, modelInput.rnCHF, modelInput.c1*T); % Reference Price
    for i=1:length(errs)
       val = PROJ_European( order,Nvec(i),alpha,r,q,T,S_0,W ,call, modelInput.rnCHF, modelInput.c1*T);
       errs(i) = log10(abs(val-ref));
    end
    
    % Plot
    plot(log2(Nvec), errs, 'r-+')
    ylabel('$log_{10}(|err|)$', 'interpreter', 'latex')
    xlabel('$log_{2}(N)$', 'interpreter', 'latex')
    title('Convergence')
    grid on;
end

if plot_smile
    figure();
    Kvec = S_0 * [0.2:.01:1.8];   % strikes to price
    
    % NOTE: there is a much more efficient version of this code for pricing many strikes
    % This script is just to provide and example of the smile
    values = zeros(1, length(Kvec));
    for i=1:length(values)
        values(i) = PROJ_European( order,N,alpha,r,q,T,S_0, Kvec(i) ,call, modelInput.rnCHF, modelInput.c1*T);
    end
    
    if call == 1
        intrinsic = max(S_0 - Kvec, 0);
    else
        intrinsic = max(Kvec - S_0, 0);
    end
    % Plot
    plot(Kvec, values)
    hold on;
    plot(Kvec, intrinsic, 'r--')
    ylabel('price')
    xlabel('strike')
    grid on;
end

