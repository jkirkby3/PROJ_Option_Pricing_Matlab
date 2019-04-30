%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MONTE CARLO ASIAN OPTION PRICER for Diffusions + Jump Diffusions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Asian options in Diffusions + Jump Diffusions
%              using the Monte Carlo Simulation
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
M    = 252;  % number of monitoring points, e.g. 252 for "daily" monitoring
Kvec = S_0*[.85 .90 .95 1 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.5 1.6];   % strikes to price

% ---------------------
% Model Params
% ---------------------
sigma = 0.2;  % diffusion parameter
jumpModel = 0;  % determines jump model, select params below

% ---------------------
% Sim Params
% ---------------------
N_sim = 10^4;  % number of simulated paths
mult = 2; % multiplier for simulation (see below) to reduce bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jumpParams = {};

if jumpModel == 1 %Normal Jumps, e.g. Merton
    lambda = 1;  muJ = -.10;  sigJ = 0.3;
    
    jumpParams.kappa = exp(muJ + .5*sigJ^2)-1;  jumpParams.lambda = lambda; jumpParams.muJ = muJ; jumpParams.sigJ = sigJ;

elseif jumpModel == 2 %Double Exponenial Jumps     
    lambda = 1;
    p_up   = 0.5; % up jump probability    
    eta1   = 25;
    eta2   = 30;
    
    kappa  = p_up*eta1/(eta1-1)+(1-p_up)*eta2/(eta2+1)-1;
    jumpParams.lambda = lambda; jumpParams.kappa = kappa; jumpParams.eta1 = eta1; jumpParams.eta2 = eta2; jumpParams.p_up = p_up;    

elseif jumpModel == 3 %Mixed normal Jumps
    lambda = 1; 
    a1 = -0.05; 
    b1 = 0.07;    
    a2 = 0.02; 
    b2 = 0.03;
    p_up = 0.6;

    kappa = p_up*exp(a1 + .5*b1^2)+ (1-p_up)*exp(a2 + .5*b2^2)  -1;
    jumpParams.lambda = lambda; jumpParams.kappa = kappa; jumpParams.a1 = a1; jumpParams.b1 = b1; jumpParams.a2 = a2; jumpParams.b2 = b2; jumpParams.p_up = p_up;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_mult = M*mult;  %time stepping to reduce bias
Spath = Simulate_Jump_Diffusion_func( N_sim, M_mult + 1, T, S_0, r, q, sigma, jumpModel, jumpParams);

disc = exp(-r*T);
[prices, stdErrs] = Price_MC_Asian_Strikes_func(Spath, call, disc, Kvec, M, mult)

% Plot
plot(Kvec, prices)
ylabel('price')
xlabel('strike')
grid on;

