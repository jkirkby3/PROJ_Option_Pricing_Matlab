function [price, opt, L] = PROJ_GMDB_DCA(proj_params, S_0, gmdb_params, r, q, modelInput)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for DCA-Style Garuanteed Minimum Withdraw Benefit (GMWB) using PROJ method
%
% Terminal Payoff:  Payoff(tau) = L*exp(g*tau) + (Gam(tau) - L*exp(g*tau))^+
%                      Gam(tau) = S_M * sum_{m=0}^M(alpha*gamma / S_m)
%                          tau  = time of death (discrete periods)
%
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
%
% NOTE: this is the SLOW "Direct" version for testing purpose. In general, use PROJ_GMDB_DCA_Fast
%
% Author: Justin Lars Kirkby
% References: 1) Equity-Linked  Guaranteed Minimum Death Benefits with Dollar Cost Averaging, J.L.Kirkby & D.Nguyen, 2021
%
% ----------------------
% Contract/Model Params 
% ----------------------
% S_0 = initial stock price (e.g. 100)
% r   = interest rate (e.g. 0.05)
% q   = dividend yield (e.g. 0.05)
% M   = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
% gmdb_params = container of GMDB contract params, see below
% modelInput =  model inputs, see below
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% proj_params = numerical params
%   proj_params.N = number of basis elements, e.g. N = 2^10
%   proj_params.L1 = gridwidth param, e.g. L1 = 8
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------
% GMDB Contract Params
% ------------------
L = gmdb_params.L;   % Guarantee Level: Set L = -1 to use ATMF value for L
alpha = gmdb_params.alpha;  % Period premium payment, paid every dt time units
gamma = gmdb_params.gamma;  % Proportion of investment retained by policyholder (fee is 1-gamma)
contract_type = gmdb_params.contract_type;  % Contract type: 1 = GMDB, 2 = GMDB-RS (Ratchet strike)
p = gmdb_params.death_prob;  % death probability distribution, must be consistent with dt
g = gmdb_params.g;

% ------------------
% Model Inputs
% ------------------
dt = modelInput.dt;  % Time increment, premiums paid / underlying S is monitored every dt
phiR = modelInput.rnCHF;  % Risk neutral CHF for time period dt


Z = gen_func(-r, dt, p);
if g == 0
    Zrg = Z;
else
    Zrg = gen_func(-(r-g), dt, p);
end

if L == -1
    MF = gen_func(r - q - g, dt, p);
    Zg = gen_func(-g, dt, p);
    L = alpha * gamma * (exp((r-q)*dt)*MF - Zg) / (exp((r-q)*dt) - 1);
end

call = 1;
N = proj_params.N;

s = 0;

for n = 1 : length(p)

    M = n;  
    T = n*dt;
    if contract_type == 2  % GMDB-RS
       W = S_0;
    else
       W = S_0*L*exp(g*T) / (alpha*gamma*(M+1));
    end
    
    pr_alpha = getTruncationAlpha(T, proj_params.L1, modelInput, proj_params.model);

    if n == 1
        W = 2*W - S_0;
        opt_v = 0.5*PROJ_European(3, N, 2*pr_alpha, r, q, T, S_0, W, call, phiR, modelInput.c1);
    else
        if contract_type == 1 || contract_type == 2  % GMDB / GMDB-RS
            ER = 0;
            opt_v = PROJ_Asian(N, pr_alpha, S_0, M, W, call, T, r, q, phiR, ER);
        elseif contract_type == 3  % European (for upper bound)
            opt_v = PROJ_European(3, N, 2*pr_alpha, r, q, T, S_0, W, call, modelInput.rnCHF_T, modelInput.c1);
            if r ~= 0  % else there is no multiplier)
                opt_v = opt_v * ((1 - exp(-r*(n+1)*dt)) / (1 - exp(-r*dt)))/(n+1);
            end
        elseif contract_type == 4  % geometric Asian (for lower bound)
            opt_v = PROJ_Geometric_Asian(N, pr_alpha, S_0, M, W, call, T, r, q, modelInput.rnSYMB);
        end
    end

    % fprintf('%.12f \n', opt_v);
    s = s + p(n) * (n + 1) * opt_v; 
    
end

opt = s * alpha * gamma / S_0;
price = L*Zrg - alpha * (exp(r*dt) - Z) / (exp(r*dt) - 1) + opt;

end


