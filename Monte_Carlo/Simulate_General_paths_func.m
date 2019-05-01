function Spath = Simulate_General_paths_func( N_sim, M, mult, T, S_0, r, q, model, modelParams, jumpModel, jumpParams )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M_mult = M*mult;  %Total number of Sampled Stock prices per path
modelType = modelParams.modelType;

if modelType == 1 % Stochastic Volatility (with Jumps)
    Spath = Simulate_StochVol_Jumps_func( N_sim, M_mult+1, T, S_0, r, q, model, modelParams, jumpModel, jumpParams);         
    
elseif modelType == 2 % Jump Diffusion (and regular diffusion, e.g. Black Scholes)
    sigma = modelParams.sigma;
    Spath = Simulate_Jump_Diffusion_func( N_sim, M_mult+1, T, S_0, r, q, sigma, jumpModel, jumpParams);    
    
elseif modelType == 3  % Stochastic Local Volatility (SLV)
    Spath = Simulate_SLV_func( N_sim, M_mult+1, T, S_0, r, q, model, modelParams);
end

end

