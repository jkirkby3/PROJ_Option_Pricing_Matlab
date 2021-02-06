function  prices = SABR_EurBarAmer_func(call, M, T, S0, Kvec, r, CTMCParams, ModParams, contract_type, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European, American, and Barrier Options using
% double Layer CTMC approximation for SABR
%
% Models Supported: SABR
% Returns: price of contract (for vector of strikes)
% Author: Justin Lars Kirkby
%
% References:  (1) General Valuation Framework for SABR and Stochastic Local Volatility
%                   Models. SIAM J. Financial Mathematics, 2018. (w/ Z. Cui
%                   and D. Nguyen)
%
% ----------------------
% Contract/Model Params 
% ----------------------
% call = 1 for call option, else Put
% contract_type = type of contact: % 1 = European, 2 = American, 3 = Down and Out Barrier
% Kvec  = strike vector
% S0 = initinal underlying value
% r   = interest rate (e.g. 0.05)
% T   = time remaining until maturity (in years, e.g. T=1)
% ModParams = model parameters: .v0, .alpha, .beta, .rho
% L =  For barrier contract, this is the barrier
%
% ----------------------
% Numerical (CTMC) Params 
% ----------------------
% CTMCParams: .m_0 = grid/state size for variance process
%             .N = grid/state stize for underlying
%             .gridMult_v = grid non-uniformity multiplier (for variance)
%             .gridMult_s = grid non-uniformity multiplier (for underlying)
%             .gamma = Grid width param for variance grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_0        = CTMCParams.m_0;
N          = CTMCParams.N;
gridMult_v = CTMCParams.gridMult_v;
gridMult_s = CTMCParams.gridMult_s;  %Grid mult param for S 
gamma      = CTMCParams.gamma;          %Grid width param for variance grid

gridMethod_v = 5;   %%% ALWAYS use 5 for this one (puts v0 on grid)
gridMethod_s = 4;   %%% 5 puts S_0 on grid, but 4 seems better (requires interpolation)

v0    = ModParams.v0;
alpha = ModParams.alpha;
beta  = ModParams.beta;
rho   = ModParams.rho;

%%%%%%%%%%%%%%%%%%%%%%
%%%   Set Asset Grid bounds
%%%%%%%%%%%%%%%%%%%%%%
if S0 < 0.5
    ls = .01*S0;   %lower bound in asset Grid (S_t) space
else
    ls = 0.001*S0;
end

us = max(4.5*S0, S0 + 10*(v0*(S0)^beta)*sqrt(T));  %upper bound in asset Grid (S_t) space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Step 1: Variance Grid / Generators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = T/M;
t = sqrt(T)/2;    %NOTE: this is different than we used to use... 

mu_func = @(u) 0*u;
sig_func = @(u) alpha*u;
mu_H = v0;
sig2_H = v0^2*(exp(alpha^2*t) - 1); 

lx = max(0.0001,mu_H - gamma*sqrt(sig2_H));
ux = mu_H + gamma*sqrt(sig2_H);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Step 1: Variance Grid / Generators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
center_v = v0;
[Q,v] = General_Q_Matrix_Newest(m_0,mu_func,sig_func,lx,ux,gridMethod_v,center_v, gridMult_v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Step 2: Asset Grid / Grid For Xtilde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = @(s) (s ).^(1-beta)/(1-beta); 
invOneBet = 1/(1-beta);

center_s = S0;   %center of asset grid (e.g. center points around the strike)
manualPoint_s = center_s;  %manually places S0 on grid
Xgrid = g(getNonUniformGrid(N, ls, us, gridMethod_s, center_s, manualPoint_s, gridMult_s)) - rho/alpha*v0;  %Grid for Xtilde

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Step 3: Generators (for Xtilde process)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nm = N*m_0;
G = zeros(Nm, Nm);  %Big Generator matrix
I = eye(N,N);

Payoff = zeros(Nm, 1);  %Terminal payoff
sqrtRho = sqrt(1-rho^2);

for j = 1:m_0 %loop through rows of big G matrix
    %%%%%%%%%%
    % Step(1): Find G_j (generator with v(j) fixed)
    %%%%%%%%%%
    nu_j = v(j);   
    muX_func_nu  = @(x) -.5*beta*(nu_j)^2.*((1-beta)*(x + rho*nu_j/alpha)).^(-1) ;  %Drift function of Xtilde with v(j) fixed
    sigX_func_nu = @(x) sqrtRho*nu_j*[x>-100];  %Constant function for each fixed variance state

    Gnu = getGenerator_Q_MatrixOnly(Xgrid, muX_func_nu, sigX_func_nu, gridMethod_v);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% FORCE absorbing vs reflecting
    Gnu(1,1) = 0; Gnu(1,2) =0;
    %Gnu(N,N) = 0; Gnu(N,N-1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%
    % Step(2): Populate the Generator matrix (recall it is block tridiagonal,
    %%%%%%%%%%
    for k = max(1,j-1):min(m_0,j+1)  %%% NOTE: we skip the ones that are known to be zeros (matrix is block tridiagonal)
        lamjk = Q(j,k);
        if j==k  %diagonal block element of G matrix
            G((j-1)*N + 1:j*N, (k-1)*N + 1:k*N ) = Gnu + lamjk*I; 
        else
            G((j-1)*N + 1:j*N, (k-1)*N + 1:k*N ) = lamjk*I;   
        end
    end 
end

%%%% Find bracketing variance gridpoint
k_0 = 2;  
while v0 >= v(k_0) && k_0 < m_0
    k_0 = k_0+1;
end
k_0 = k_0 - 1;  %left bracketing point:  v(k_0) <= v0 < v(k_0 +1)

%%%% Find bracketing Xtilde gridpoint
x0 = g(S0) - rho*v0/alpha;  %initial value on Xtilde grid for (S_0, v0)
j_0 = 2;
while x0 >= Xgrid(j_0) && j_0 < N
    j_0 = j_0+1;
end
j_0 = j_0 - 1;  %left bracketing point:   Xtilde(j_0) <= x0 < Xtilde(j_0 +1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VALUE : using recursive method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = expm(G*dt);   %transition matrix of one dimensional CTMC
initialIndex = (k_0-1)*N + j_0;  %index corresponding to initial conditions

prices = zeros(length(Kvec),1);

for k = 1:length(Kvec)
    K = Kvec(k);
    %%% Calculate Payoff
    if call == 1 
        for j = 1:m_0
            Payoff((j-1)*N +1:j*N) = max(0, ((1-beta)*(max(0,Xgrid + rho*v(j)/alpha) )).^invOneBet - K) ;  %Portion of payoff corresponding to v(j)
        end
    else 
       for j = 1:m_0
            Payoff((j-1)*N +1:j*N) = max(0, K - ((1-beta)*(max(0,Xgrid + rho*v(j)/alpha) )).^invOneBet );  %Portion of payoff corresponding to v(j)
       end
    end
    %%% Now Price
    if contract_type == 3  % Down and out for now
        %determine which states remain alive
        alive = zeros(Nm,1);
        for j = 1:m_0
            cons = g(L) - v(j)*rho/alpha;
            alive((j-1)*N +1:j*N) = (Xgrid > cons);
        end
        pVec = alive.*Payoff;
        for m = M-1:-1:0
            pVec = exp(-r*dt)*alive.*(P*pVec);     
        end
    else %% either American or European
        pVec = exp(-r*dt)*(P*Payoff);  %initialize the value (at last period)
        if contract_type == 2
            for m = M-2:-1:0
                pVec = exp(-r*dt)*(P*pVec);  %continuation value
                pVec = max(pVec, Payoff);   %max of continuation and intrinsic value
            end
        elseif contract_type == 1
            for m = M-2:-1:0
                pVec = exp(-r*dt)*(P*pVec);  %continuation value
            end
        end
    end
    
    if gridMethod_s == 5  %Both v0 and x0 are on grid  (we assume v0 is a member of vol grid, ie as long as gridMethod_v = 5)
        prices(k) = pVec(initialIndex);    
    elseif gridMethod_s == 4 %x0 is not on grid (though we assume v0 is, ie as long as gridMethod_v = 5)
        price1 = pVec(initialIndex);    %corresponds to j_0  
        price2 = pVec(initialIndex+1);  %corresponds to j_0 + 1
        prices(k) = price1 + (price2 - price1)*(x0 - Xgrid(j_0))/(Xgrid(j_0+1) - Xgrid(j_0));  %LINEAR INTERPOLATION
    end
end
    
end

