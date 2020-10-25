function [vals, c_index_1, c_index_2, y_1, y_2] = price_2d_ctmc( S_0s, T, r, rho, sigmas, qs, params, contractParams, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European/Bermudan/Barrier Options using CTMC approximation method
% Models Supported: 2D Diffusions
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% References:  (1) A General Continuous Time Markov Chain Approximation for
%               Multi-Asset option pricing with systems of correlated diffusions,
%               Applied Math. and Comput., 2020 (JL Kirkby, Duy Nguyen, Dang Nguyen)
%
% ----------------------
% Model Params 
% ----------------------
% S_0s   = initial asset prices
% T      = time remaining until maturity (in years, e.g. T=1)
% r      = interest rate (e.g. 0.05)
% rho    = correlation between brownian motions 
% sigmas = volatilities per asset
% qs     = div yeilds per asset
% M      = num monitoring points for Barrier/Bermudan (also controls num steps for multi-step European pricing)
%
% ----------------------
% Contract Params  (contractParams)
% ----------------------
%
% contractParams.payoff_type:
% 
% 1: Linear, G = S_1  (linear payoff in first underlying)
% 2: Linear, G = S_2  (linear payoff in second underlying)
% 3: Exchange, G = (S_1 - S_2)^+
% 4: Spread,  G = (S_1 - S_2 - K)^+   (NOTE: must set strike, K)
% 5: Geometric Basket Call / Put,  G = (sqrt(S_1) * sqrt(S_2) - K)^+  (for the call)
% 6: Arithmetic Basket Call / Put,  G = (sqrt(S_1) * sqrt(S_2) - K)^+  (for the call)
% 7: Call-on-Max and Put-on-Min, Gcall = (max(S_1,S_2) - K)^+ , Gput = (K - min(S_1,S_2))^+
% 8: Call/put on just S_2, G = (S_2 - K)^+  (for the call)
% 9: Best-of / Worst-of,  G = max(S_1,S_2), G = min(S_1,S_2)
% 
% contractParams.contract:
%
% 1: European, Single Step Pricing
% 2: European, Multi Step Pricing (M above controls number of steps)
% 3: Bermudan (M above controls number of monitoring points)
% 4: Barrier (M above controls number of monitoring points)
%
% Note: for barrier option, set:
%       contractParams.barriers_1 = lower/upper barriers on first asset  (e.g. [0 50])
%       contractParams.barriers_2 = lower/upper barriers on second asset  (e.g. [0 50000000000] to disable barrier on S_2)
%
% ----------------------
% Numerical (CTMC) Params 
% ----------------------
% params = CTMC parameters
% params.m_0 = num CTCM states
% params.num_devs = num std devs used in the grid
% params.gridMethod = choose the grid method (several RnD versions)
% params.GridMultParam = non-uniformity param, in (0,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 9
    M = 1; % M is only needed for 
end
dt = T/M;

contract = contractParams.contract;

if contract == 1  % European
    dt = 1; M = 1;
end

method = 4;
num_devs = params.num_devs;
m_0 = params.m_0;
GridMultParam = params.GridMultParam;
gridMethod = params.gridMethod;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drifts = r - qs;

R = [1 rho; rho 1];  % Correlation Matrix
[ L, D, C, Cinv ] = get_transform_matrices_2d( R, method );

% Now Define New Uncorrelated System  (the dc underscore)
[drift_dc, sigma_dc ] = decorrelate(sigmas, drifts, C, D );

[Ls_dc, Rs_dc ] = get_CTMC_decorr_boundaries(sigmas, C, T, 0, sigma_dc, num_devs);
Y_0s = [0 0];
    
% Form CTMC 1
center = Y_0s(1);
mu_func = @(s) drift_dc(1)*[s>-100000];
sig_func = @(s) sigma_dc(1)*[s>-100000];
[Q, y_1, c_index_1] = Q_Matrix(m_0, mu_func,sig_func,Ls_dc(1),Rs_dc(1),gridMethod,center, GridMultParam);
P1 = expm(Q*dt);

% Form CTMC 2
center = Y_0s(2);
mu_func = @(s) drift_dc(2)*[s>-100000];
sig_func = @(s) sigma_dc(2)*[s>-100000];
[Q, y_2, c_index_2] = Q_Matrix(m_0, mu_func,sig_func,Ls_dc(2),Rs_dc(2),gridMethod,center, GridMultParam);
P2 = expm(Q*dt);

G = get_payoff_G_matrix_from_ygrid_2d( y_1, y_2, S_0s, sigmas, rho, contractParams);

if contract == 1  % European, Price by single step
    vals = exp(-r*T)*P1*G*P2.';
    
elseif contract == 2  % European, Price By Multi Step
    vals = G;
    for m=M-1:-1:0
        vals = exp(-r*dt)*P1*vals*P2.';
    end 
    
elseif contract == 3  % Bermudan
    vals = G;
    for m=M-1:-1:0
        vals = max(exp(-r*dt)*P1*vals*P2.', G);
    end
    
elseif contract == 4 % Barrier
    b1 = contractParams.barriers_1; L1 = b1(1); U1 = b1(2);
    b2 = contractParams.barriers_2; L2 = b2(1); U2 = b2(2);
    B = ones(m_0, m_0);
    for i = 1:m_0
        y1 = y_1(i);
        S1 = S_0s(1)*exp(sigmas(1)*y1);
        if S1 < L1 || S1 > U1
            B(i,:) = 0;
        else
            for j = 1:m_0
                y2 = y_2(j);
                S2 = S_0s(2)*exp(sigmas(2)*(y2 + rho*y1));
                if S2 < L2 || S2 > U2
                    B(i,j) = 0;
                end
            end
        end
    end
    vals = G.*B;
    for m=M-1:-1:0
        vals = exp(-r*dt)*(B.*P1*vals*P2.');
    end 
end

end

