%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price European options under Regime Switching Diffusion Models using Monte Carlo simulation
% Author:      Justin Kirkby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);

addpath('../')

% ---------------------
%  Contract/Market Params
% ---------------------
call = 1;    %For call use 1 (else, its a put)
S_0  = 100;  %Initial price
r    = .05;  %Interest rate
q    = .00;  %dividend yield
T    = 1;    %Time (in years)
Kvec = S_0*[.85 .90 .95 1 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.5 1.6];   % strikes to price

% ---------------------
% Regime Switching Diffusion Params
% ---------------------
% Transition Matrix (dictates how the regimes transition)
Q = [-1 0.5 0.5;
    0.5 -1 0.5; 
    0.5 0.5 -1];  

drift_vec = [r-q  r-q  r-q];  % Drift in each state
sigma_vec = [0.15  0.25  0.35]; % Volatility in each state

initial_state = 1;

% ---------------------
% Sim Params
% ---------------------
N_sim = 10^5;
M = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Spath = Simulate_RegimeSwitching_Diffusion_func( N_sim, M, T, S_0, drift_vec, sigma_vec, Q, initial_state);
histogram(Spath(:,end))

disc = exp(-r*T);
[prices, stdErrs] = Price_MC_European_Strikes_func(Spath, disc, call, Kvec )

plot(Kvec, prices)
ylabel('price')
xlabel('strike')
grid on;


