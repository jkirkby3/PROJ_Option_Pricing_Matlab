%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ASIAN OPTION PRICER  (RUN SCRIPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Gauranteed Minimum Death Benefits (GMDB) in Levy Models using the PROJ method
%               This version is based on a dollar cost average style investment account (see reference paper below)
%
% Terminal Payoff:  Payoff(tau) = L*exp(g*tau) + (Gam(tau) - L*exp(g*tau))^+
%                      Gam(tau) = S_M * sum_{m=0}^M(alpha*gamma / S_m)
%                          tau  = time of death (discrete periods)
%                            M  = number of periods until time of death (each period length dt)
%
% Author:      Justin Kirkby
% Reference:    1) Equity-Linked Guaranteed Minimum Death Benefits with Dollar Cost Averaging, J.L.Kirkby & D.Nguyen, 2021
%               2) An Efficient Transform Method For Asian Option Pricing, SIAM J. Financial Math., 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);
addpath('../RN_CHF')
addpath('../Helper_Functions')
addpath('../Asian_Options')
addpath('../European_Options')
addpath('../Geometric_Asian_Options')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 1) CHOOSE CONTRACT/GENERAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_0  = 100;  %Initial price
r    = 0.01;  %Interest rate
q    = 0.00;  %dividend yield
age  = 30;  % age at time of purchase/valuation

gmdb_params = {};
gmdb_params.contract_type = 1;  % 1 = GMDB, 2 = GMDB-RS (ratchet strike)
gmdb_params.alpha = 2*S_0;      % Periodic investment amount
gmdb_params.gamma = 0.92;       % gamma in [0,1], Administrative fee is (1-gamma)
gmdb_params.L = 10000;          % Gaurantee Amount at death, set to -1 to take ATMF value
gmdb_params.g = 0.01;           % Growth rate of guarantee, becomes L*exp(g*tau)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Step 2) CHOOSE MODEL PARAMETERS (Levy Models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------
% Choose Model for death/mortality
% -------------
death_model = 1; % 1 = Mortality Table, 2 = Combination of exponentials (see reference paper above)

if death_model == 1
    death_prob = make_mortality_table_pmf( age );
elseif death_model ==2
    death_prob = make_combo_2_expos_pmf( 3, -2, 0.08, 0.12, 110 - age + 1);
elseif death_model ==3 
    % manually specify probability of death
    death_prob = [0.01 0.05 0.1 0.2 0.1 0.05 0.001 0.5 0.2];
    death_prob = death_prob / sum(death_prob);
end
gmdb_params.death_prob = death_prob;

% -------------
% Choose Levy Model for risky asset
% -------------
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
freq = 1;   % Harcoded: frequency of payments (must correspond to death probability mass function)
dt = 1/freq;

proj_params = {};
proj_params.N = 2^7;  % number of PROJ grid points
proj_params.L1 = 8;  % determines grid width (usually set L1 = 8 to 15 for Levy)
proj_params.model = model;

modelInput = getModelInput(model, dt, r, q, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_death_prob = 1; 
if plot_death_prob == 1
    plot(1:length(death_prob), death_prob, 'k-', 'linewidth', 1.1)
    xlabel('$M_\tau = n$', 'interpreter', 'latex')
    ylabel('probability, $p^\omega_n$', 'interpreter', 'latex')
end


tic
proj_params.L1 = 12;  % determines grid width (usually set L1 = 8 to 15 for Levy)
[price2, opt2, L2] = PROJ_GMDB_DCA_Fast(proj_params, S_0, gmdb_params, r, q, modelInput);
toc

fprintf('Price: %.8f, Option: %.8f, L: %.8f \n', price2, opt2, L2)
