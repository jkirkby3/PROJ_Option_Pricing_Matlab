function price = TC_Levy_VarianceSwap(r,q, T, M, levyExponent, hFunc, n, ParamsDiffus, ParamsCtmc)
% levyExponent is function handle: e.g. -1/2*i*xi - 1/2*xi^2 for case of Heston as time changed levy process
% hFunc is function handle: tau = int_0^T h(X_s)ds
% n is number of discrete points for the discrete version, set n = 0 for continuous
% ParamsDiffus are the diffusion parameters of X_s, where tau = int_0^T h(X_s)ds

% Martingale Adjustment for exponent: ensure the driving Levy process itself is a martingale
% Note, we model the RN drift (r-q)*T separately... then price using PROJ 
levyExponentRN = @(u)levyExponent(u) - 1i*levyExponent(-1i)*u;  

dt = T/M;

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
boundaryMethod = 1;
[G,stateGrid] = Q_Matrix_AllForms(Nx,mu_func,sig_func,lx,ux,gridMethod, varGridMult, center, boundaryMethod);

H = hFunc(stateGrid);  % This will be used along diagonal

%%%%////////////////////////////////////////////////////////
%%%% Find bracketing states
%%%%////////////////////////////////////////////////////////
j_0 = 2;  %k_0 will satisfy v(k_0) <= v0 < v(k_0 +1)
while v0 >= stateGrid(j_0) && j_0 < Nx
    j_0 = j_0+1;
end
j_0 = j_0 - 1;

%%%%////////////////////////////////////////////////////////
% Risk Neutral CHF of Log Return: ln(S_T/S_0)
%%%%////////////////////////////////////////////////////////
RNdrift = (r-q)*dt;
delt = 5e-03;
z = [-delt delt];   % To evaluate derivative of chf

if n == 0  % Continuous case
    chf = rnCHF_cont(z, levyExponentRN, H, G, dt, stateGrid, RNdrift);
else   % Discrete case
    chf = rnCHF_disc(z, levyExponentRN, H, G, dt, n, stateGrid, RNdrift);
end

Evec = zeros(Nx, 1);

for j = 1:Nx
    M_mh = chf(j,1);
    M_0 = 1;
    M_ph = chf(j,2);
    Evec(j) = -(M_ph - 2*M_0 + M_mh) /delt^2;
end

Pm = expm(dt*G);
PmPm = Pm;
sumPm = eye(Nx, Nx);
for m = 2:M
    sumPm = sumPm + PmPm;
    if m < M
        PmPm = Pm*PmPm;
    end
end

prices = sumPm*Evec; 
price = prices(j_0)/T;

end

