%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Variance Swaps - Time Changed Heston Option Pricer (This example is a time-changed representation of Heston's model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to Price Discrete Variance Swap under Hestons Model
% Author:      Justin Kirkby
% References:  (1) A General Framework for tim changed Markov Processes and Applications
%              European J. Operational research, 2019
%              (2) Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform, SIAM J. Financial Math., 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[folder, name, ext] = fileparts(which( mfilename('fullpath')));
cd(folder);

addpath('../Helper_Functions')
addpath('../../STOCHASTIC_VOL/DV_Swaps_Options/Analytical_Swaps')
addpath('../../STOCHASTIC_VOL/Helper_Functions')

% ==========================
% CTMC Params
% ==========================
ParamsCtmc.varGridMult = .05;
ParamsCtmc.gamma = 6;  % Heston gamma = 4 is good for T ~ 1
ParamsCtmc.Nx = 70; %the number of Markov states
n = 0; % number of time steps in time disretization... set to 0 to do continuous time version

% ==========================
% Contract/Market Params
T = 1;%time to maturity
r = 0.05;  % interest rate
q = .00;   % dividend yeild
% ==========================

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Heston Example   % Y_t = W_t - t/2   
eta = 3;  
theta = 0.04;
Sigmav = 0.1;
v0 = 0.04;
rho = 0.0;  % NOTE: rho must be zero, time-changed approach here requires zero correlation (otw its approximation)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParamsDiffus.model = 1;
ParamsDiffus.eta = eta; ParamsDiffus.rho = rho; ParamsDiffus.theta = theta; ParamsDiffus.Sigmav = Sigmav; ParamsDiffus.v0 = v0;

hFunc = @(u) u;   % tau = int h(X_s) ds
levyExponent = @(z) -0.5*1i*z - 0.5*z.^2; 

mvec = [12 52 252];  % differnt values of M, the number of discete monitoring dates
time_iter = 20;

for m = 1:length(mvec)
    M = mvec(m);
    tic
    for i = 1:time_iter
        price = TC_Levy_VarianceSwap(r, q, T, M, levyExponent, hFunc, n, ParamsDiffus, ParamsCtmc);
    end
    time = toc/time_iter;
    
    % ANALYTICAL
    [ref, KcH] = hestonfairstrike(r, v0, theta, eta, Sigmav, T, rho, M);
    error = price - ref;
    
    fprintf('%.0f & %.8f & %.8f & %.1e & %.3f \n', m, ref, price, abs(error), time)
    
end



