%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GMDB OPTION PRICER with Surrender
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Script to Price Gauranteed Minimum death benefits with period fees and early surrender in Levy Models
%              using the PROJ method
% Author:      Justin Kirkby
% References:  (1) Valuation and optimal surrender of variable annuities
%                   with guaranteed minimum benefits and periodic fees, 
%                   Kirkby and Aguilar 2022, Scandinavian Actuarial Journal

%              (2) Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform, SIAM J. Financial Math., 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../')
addpath('../../RN_CHF')
addpath('../../Helper_Functions')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r    = .02;  % interest rate
q    = .01;  % dividend yield

age = 30;      % age of policyholder at time of purchase

gmdb_params = {};
gmdb_params.F_0 = 1;            % initial fund value
gmdb_params.alpha_fee = 0.02;   % period fee rate
gmdb_params.gamma = 1.0;        % surrender penalty rate, 1.0 = 100% fund lost upon surrender, 0.0 = 0% lost, ie no penalty
gmdb_params.g = 0.03;           % floor on growth
gmdb_params.c = 0.3;            % cap on growth

Ts = [5, 10, 20, 30, 50];  % Years until maturity (will price each of these)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS  (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    params.alpha = 9.5;
    params.beta  = -0.4;
    params.delta = 2;
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 3) CHOOSE PROJ PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------
% APPROACH 1: Cumulant Based approach for grid width
% (see "Robust Option Pricing with Characteritics Functions and the BSpline Order of Density Projection")
%---------------------
logN  = 12;   %Uses N = 2^logN  gridpoint 
L1 = 11;  % determines grid witdth (usually set L1 = 8 to 15 for Levy)

N = 2^logN;    % grid roughly centered on [c1 - alph, c1 + alph]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

price_refs = zeros(1, length(Ts));
price_projs = zeros(1, length(Ts));
times = zeros(1,length(Ts));
errs = zeros(1, length(Ts));
reps = 5;

for i=1:length(Ts)
    T = Ts(i);
    M = T;
    %%%  Note: rnCHF is the risk netural CHF, c1,c2,c4 are the cumulants
    modelInput = getModelInput(model, T/M, r, q, params);
    params.zeta = -modelInput.rnCHF(-1i);
    
    death_prob = make_mortality_table_pmf( age );
    death_prob = death_prob(1:M);
    death_prob(M) = death_prob(M) + (1 - sum(death_prob));

    %%% Conditional death prob for dynamic programming
    death_prob_cond = make_mortality_table_pmf( age,1);
    death_prob_cond = death_prob_cond(1:M);
    death_prob_cond(M) = death_prob_cond(M) + (1 - sum(death_prob_cond));
    
    gmdb_params.death_prob = death_prob;
    gmdb_params.death_prob_cond = death_prob_cond;

    alpha = getTruncationAlpha(T, L1, modelInput, model);
    
    if model == 1 && gmdb_params.gamma == 1.
        price_refs(i) = Price_GMDB_BSM_NoSurrender(T, M, gmdb_params, params.sigmaBSM, r, q);
    else
        price_refs(i) = PROJ_GMXB_Surrender(T, M, gmdb_params, modelInput, 2*N, alpha);
    end
    
    tic
    for k = 1:reps
        price_projs(i) = PROJ_GMXB_Surrender(T, M, gmdb_params, modelInput, N, alpha);
    end
    times(i) = toc / k;
    errs(i) = abs(price_refs(i) - price_projs(i));
end

fprintf('$T$ & %.0f & %.0f & %.0f & %.0f & %.0f \\\\ \n', Ts)
fprintf('Ref. & %.8f & %.8f & %.8f & %.8f & %.8f \\\\ \n', price_refs)
fprintf('PROJ & %.8f & %.8f & %.8f & %.8f & %.8f \\\\ \n', price_projs)
fprintf('Err. & %.2e & %.2e & %.2e & %.2e & %.2e \\\\ \n', errs)
fprintf('Time(s) & %.3f & %.3f & %.3f & %.3f & %.3f \\\\ \n', times)
