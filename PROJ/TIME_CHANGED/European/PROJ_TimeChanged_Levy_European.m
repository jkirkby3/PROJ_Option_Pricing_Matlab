function price = PROJ_TimeChanged_Levy_European(r,q,S_0,T,W,call, levyExponent, hFunc, n, ParamsDiffus, ParamsCtmc, ParamsProj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European Options under time-changed Levy Models,
%        using Markov-Chain Approximation + PROJ method
% Models Supported: Time Changed Levy Processes, such as Heston's Model (with rho=0)
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% References:   1) A General Framework for Time-Changed Markov Processes and Applications. 
%                  European J. Operational Research (2018), Cui, Z., Kirkby, J.L., and Nguyen, D.
%
%               2) Efficient Option Pricing by Frame Duality with the Fast Fourier Transform. 
%                  SIAM J. Financial Math (2015), Kirkby, J.L
% ----------------------
% Contract/Model Params 
% ----------------------
% S_0 = initial stock price (e.g. 100)
% W   = strike  (e.g. 100)
% r   = interest rate (e.g. 0.05)
% q   = dividend yield (e.g. 0.05)
% T   = time remaining until maturity (in years, e.g. T=1)
% call  = 1 for call (else put)
% hFunc = function handle, the integrand in time change, tau = int_0^T h(X_s)ds
% levyExponent = function handle: e.g. -1/2*i*xi - 1/2*xi^2 for case of Heston as time changed levy process
% n = number of discrete points for the discrete version, set n = 0 for continuous
% ParamsDiffus =  diffusion parameters of X_s, where tau = int_0^T h(X_s)ds
%
% ----------------------
% Numerical (CTMC/PROJ) Params 
% ----------------------
% ParamsCtmc = container of CTMC approximation params
% ParamsProj = container of PROJ params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Martingale Adjustment for exponent: ensure the driving Levy process itself is a martingale
% Note, we model the RN drift (r-q)*T separately... then price using PROJ 
levyExponentRN = @(u)levyExponent(u) - 1i*levyExponent(-1i)*u;  


%%% Set CTMC Parameters (CTMC approximates the diffusion)
gridMethod = 6;  % NOTE: will need to interpolate below if gridMethod != 6
varGridMult = ParamsCtmc.varGridMult;
gamma = ParamsCtmc.gamma;    
Nx = ParamsCtmc.Nx; %the number of Markov states

%%% The time change driver process X_t is approximated by CTMC
if ParamsDiffus.model == 1 %CIR , e.g. HESTON  (eta, theta, Rho, Sigmav, v0)
    t = T;
    eta = ParamsDiffus.eta; theta = ParamsDiffus.theta;  Sigmav = ParamsDiffus.Sigmav; v0 = ParamsDiffus.v0;  % Rho = ParamsDiffus.rho;
    mu_func  = @(v)eta*(theta-v);  %variance process drift
    sig_func = @(v)Sigmav*sqrt(v); %variance process vol
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmav^2/eta*v0*(exp(-eta*t)-exp(-2*eta*t)) +theta*Sigmav^2/(2*eta)*(1-exp(-eta*t)+exp(-2*eta*t));
    
    lx = max(0.00001,mu_H - gamma*sqrt(sig2_H));  %variance grid lower bound
    ux = mu_H + gamma*sqrt(sig2_H);  %diffusion grid upper bound
end

center = v0; 
[G, stateGrid]  =  Q_Matrix_AllForms(Nx,mu_func,sig_func,lx,ux,gridMethod, varGridMult, center);

H = hFunc(stateGrid);  % This will be used along diagonal

%%%%////////////////////////////////////////////////////////
%%%% Find bracketing states
%%%%////////////////////////////////////////////////////////
k_0 = 2;  %k_0 will satisfy v(k_0) <= v0 < v(k_0 +1)
while v0 >= stateGrid(k_0) && k_0 < Nx
    k_0 = k_0+1;
end
k_0 = k_0 - 1;

%%%%////////////////////////////////////////////////////////
% Risk Neutral CHF of Log Return: ln(S_T/S_0)
%%%%////////////////////////////////////////////////////////
if n == 0  % Continuous case
    rnCHF = @(z)exp(1i*(r-q)*T*z).*rnCHF_cont_single(z, levyExponentRN, H, G, T, k_0, stateGrid);
else   % Discrete case
    rnCHF = @(z)exp(1i*(r-q)*T*z).*rnCHF_disc_single(z, levyExponentRN, H, G, T, n, k_0, stateGrid);  % Makes risk neutral
end

c1 = 0;
order = ParamsProj.order;
alph = ParamsProj.alph;
N_proj = ParamsProj.N_proj;

price = PROJ_European(order, N_proj, alph, r, q, T, S_0, W, call, rnCHF, c1);

end

function chf = rnCHF_cont_single(z, levyExponentRN, H, G, T, k_0, stateGrid)
    numStates = length(stateGrid);
    numZ = length(z);
    
    ones_ = ones(numStates, 1);   % column vector
    chf = zeros(numStates, numZ); % chf for each initial state 
    
    for j = 1:numZ   % populate one column (one theta(i) value) at a time
        chf(:, j) = expm( T*(G + diag(levyExponentRN(z(j))*H))) * ones_;   % to initialize matrix vector mult
    end
    
    chf = chf(k_0,:);
end

function chf = rnCHF_disc_single(z, levyExponentRN, H, G, T, n, k_0, stateGrid)
    % k_0 is the bracketing index: grid(k_0) <= v0 < grid(k_0+1)
    %Returns the chf with is linear interpolation of the two bracketing chfs
    dt = T/n;
    P = expm(dt*G);
    
    numStates = length(stateGrid);
    numZ = length(z);
    
    ones_ = ones(numStates, 1);   % column vector
    chf = zeros(numStates, numZ); % chf for each initial state 

    for j = 1:numZ   % populate one column (one theta(i) value) at a time
        E = diag(exp(levyExponentRN(z(j))*H*dt));
        chf(:, j) = E*ones_;   % to initialize matrix vector mult
        EP = E*P;
        for k = 1:n
            chf(:, j) = EP*chf(:, j); 
        end
     
    end
    
    chf = chf(k_0,:);
end



