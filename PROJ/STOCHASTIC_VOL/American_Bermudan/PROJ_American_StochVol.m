function price = PROJ_American_StochVol(numeric_param, M, r, T, S_0, W, psi_J, model, modparam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for American PUT Option using CTMC Approximation + PROJ method
% Models Supported: Stochastic Volatility (including jumps)
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% References:  (1) A unified approach to Bermudan and Barrier options under stochastic
%               volatility models with jumps. J. Economic Dynamics and Control, 2017
%              (2) Robust barrier option pricing by Frame Projection under
%               exponential Levy Dynamics. Applied Mathematical Finance, 2018.
%
% ----------------------
% Contract Params 
% ----------------------
% T  : number of years (T = 2 is two years, T = .5 is half a year)
% M  : number of subintervals of [0,T] (total of M+1 points in time grid)
% W  : strike  (used instead of K)
%
% ----------------------
% Model Params 
% ----------------------
% S_0: initial Underlying
% r  : interest rate 
% psi_J: characteristic exponenent of jump part...
%        function handdle: psi_J(xi) = lambda*(phi(xi) -1)
% model:
%        1 = HESTON:      Sigmav, v0, rho, eta, theta
%        2 = STEIN-STEIN: Sigmav, v0, rho, eta, theta
%        3 = 3/2 MODEL:   Sigmav, v0, rho, eta, theta
%        4 = 4/2 MODEL:   Sigmav, v0, rho, eta, theta, aa, bb
%        5 = HULL-WHITE:  Sigmav, v0, rho
%        6 = SCOTT:       Sigmav, v0, rho, eta, theta
%        7 = ALPHA-HYPER: Sigmav, v0, rho, eta, theta
%
% modparam: contains all necessary params for the specific model (see below during assingment which ones are needed)
%
% ----------------------
% Numerical Params 
% ----------------------
% numeric_parm: container of numerical params
%   N  : size of density grid (value grid is K:=N/2)
%   alph: density gridwith param, density on [-alph,alph]... value grid width = alph
%   m_0: number of states to approximate the Heston model with
%   gamma: var grid width parameter, grid is +/- gamma*stddev(variance process)
%   gridMethod: which type of var grid to use (typcially use 4)
%-------------------------------
%
%%% Note: be careful about the parameter Rho (vs rho used in algorithm)

N = numeric_param.N;
alph = numeric_param.alph;
m_0 = numeric_param.m_0;
gridMethod = numeric_param.gridMethod;
gamma = numeric_param.gamma;
gridMultParam = numeric_param.gridMultParam;

K    = N/2;
dx   = 2*alph/(N-1);
lws  = log(W/S_0);
dt   = T/M;

%%% GRID which aligns 0 as well as log(W/S_0)
nnot = K/2;
dxtil = dx; %for now... change to dxtil = 2*alpha/(N-1)
nbar = floor(lws/dx +K/2);
if abs(lws)<dxtil
    dx = dxtil;  
elseif lws<0
    dx = lws/(1+nbar - K/2);
    nbar = nbar+1;
elseif lws>0
    dx = lws/(nbar-K/2);
end
a = 1/dx;
xmin = (1-K/2)*dx;


%%%%////////////////////////////////////////////////////////
%%%% Initialize THETA ... 
%%%%////////////////////////////////////////////////////////
THET  = zeros(K,m_0);  %Each column corresponds to a regime state

%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;


%%%% PAYOFF CONSTANTS-----------------------------------
varthet_01 = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
varthet_m10 = exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4)/18;
varthet_star = varthet_01 + varthet_m10;

%%%%%%%
Gs  = zeros(K,1);
Gs(1:nbar) = exp(xmin + dx*(0:nbar-1))*S_0;

%Define Terminal Theta Coeffs 
ThetM = zeros(K,1);
ThetM(1:nbar-1) = W - varthet_star*Gs(1:nbar-1);
ThetM(nbar)     = W*(.5 - varthet_m10);
Gs(1:nbar)      = W - Gs(1:nbar); %For American options (otherwise no need for Gs to be defined separately)


%%%%////////////////////////////////////////////////////////
%%%% Intialize Q matrix and variance set
%%%%////////////////////////////////////////////////////////
t = T/2;
[lx, v0, ux] = get_variance_grid_boundaries( model, modparam, t, gamma);

[mu_func,  sig_func] = get_SV_variance_grid_diffusion_funcs( model,  modparam);
boundaryMethod = 1;

center = v0; %this is where grid clusters... we can experiment with other choices.. 

[Q,v]  = Q_Matrix_AllForms(m_0,mu_func,sig_func,lx,ux,gridMethod, gridMultParam, center, boundaryMethod);

%%%%////////////////////////////////////////////////////////
%%%% Populate the Matrix Exponentials
%%%%////////////////////////////////////////////////////////
dxi    = 2*pi*a/N;
xi     = dxi*(0:N-1)';  

[v1, v2, fv] = get_SV_matrix_expo_inputs( model,  modparam, psi_J, dt, v, dxi, r);
EXP_A = get_SV_matrix_exponential( Q, dt, xi, v1, v2, fv, psi_J, m_0, N );

%%%%////////////////////////////////////////////////////////
%%%% Construct Toepliz Array Of Arrays
%%%%////////////////////////////////////////////////////////
a2    = a^2;  
Cons2 = 24*a2*exp(-r*dt)/N;
zmin  = (1 - K)*dx;  %Kbar corresponds to zero

xi     = dxi*(1:N-1); %REDFINED FROM ABOVE
hvec = exp(-1i*zmin*xi).*(sin(xi/(2*a))./xi).^2./(2+cos(xi/a));

BETA = zeros(N,m_0,m_0);  %to access the (j,k)th toeplitz data, use BETA(:,j,k)
grand = zeros(1,N-1);
      
%%% NOTE the (k,j) rather than (j,k)
for j=1:m_0
    for k = 1:m_0
        for n=1:N-1  % MAKE this a .* if possible
            grand(n) = hvec(n)* EXP_A(k,j,n+1);  %n+1 since EXP_A has the first element defined
        end
        beta  = Cons2*real(fft([EXP_A(k,j,1)/(24*a2) grand]));   %%%%  NOTE: all toep matrices incorporate exp(-r*dt)
        toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];   
        BETA(:,j,k) = fft(toepM);
    end
end

%%%%////////////////////////////////////////////////////////
%%%% Initialize Continuation Value
%%%%////////////////////////////////////////////////////////
% CONT = repmat(p(1:K),1,m_0);  %continuation value in each state (initialize with all equal)
CONT = zeros(K,m_0);
ThetTemp = fft([ThetM(1:K);zeros(K,1)]);
for j=1:m_0
    for k = 1:m_0
        p = ifft(BETA(:,j,k).*ThetTemp);  
        CONT(:,j) = CONT(:,j)+ p(1:K);
    end
end

%%%%////////////////////////////////////////////////////////
%%%% LOOP through time
%%%%////////////////////////////////////////////////////////

kstr_vecInit = nbar*ones(m_0,1);   %in general one per state

for m = M-2:-1:0
    kstr_vec = kstr_vecInit;  %%%% FOR NOW!!! RESET EACH TIME (in future, we can start from previous known value to save cost)
    
    %Step 1: update THETA
    for j=1:m_0
        while kstr_vec(j)>2 && CONT(kstr_vec(j),j)> Gs(kstr_vec(j))
            kstr_vec(j) = kstr_vec(j) -1;
        end
        if kstr_vec(j)>=2
            xkstr = xmin +(kstr_vec(j) -1)*dx;

            Ck1 = CONT(kstr_vec(j)-1,j); Ck2 = CONT(kstr_vec(j),j); Ck3 = CONT(kstr_vec(j)+1,j);
            %%% Linear interp of payoff
            Gk2 = Gs(kstr_vec(j)); Gk3 = Gs(kstr_vec(j)+1);
            %%% xstr
            tmp1 = Ck2-Gk2; tmp2 = Ck3 - Gk3;
            xstrs  = ((xkstr+dx)*tmp1 - xkstr*tmp2)/(tmp1-tmp2);
        else
            kstr_vec(j) = 1; xstrs = xmin; xkstr = xmin; 
        end
           
        rho = xstrs-xkstr;
        zeta = a*rho;
        
        zeta2 = zeta^2; zeta3 = zeta*zeta2; zeta4 = zeta*zeta3;
        zeta_plus = zeta*q_plus; zeta_minus = zeta*q_minus;
        rho_plus = rho*q_plus; rho_minus = rho*q_minus;

        ed1 = exp(rho_minus); ed2 = exp(rho/2); ed3 = exp(rho_plus);
        dbar_1 = zeta2/2;
        dbar_0 = zeta - dbar_1;         %  dbar_1 = zeta + .5*((zeta - 1)^2 - 1);
        d_0    = zeta*(5*( (1-zeta_minus)*ed1 + (1-zeta_plus)*ed3 ) + 4*(2-zeta)*ed2)/18;
        d_1    =  exp(-dx)*zeta*( 5*(zeta_minus*ed1 + zeta_plus*ed3) + 4*zeta*ed2 )/18;    
        
        THET(1:kstr_vec(j) -1,j) =  ThetM(1:kstr_vec(j)-1);
        Ck4 = CONT(kstr_vec(j) +2,j);
        
        THET(kstr_vec(j),j) = W*(.5 + dbar_0) -S_0*exp(xkstr)*(varthet_m10 +d_0)...
             + zeta4/8*(Ck1 -2*Ck2 + Ck3) +zeta3/3*(Ck2-Ck1)...
             + zeta2/4*(Ck1 +2*Ck2 - Ck3) - zeta*Ck2...
             -Ck1/24 +5/12*Ck2 +Ck3/8;
     
        THET(kstr_vec(j)+1,j) = W*dbar_1 - S_0*exp(xkstr+dx)*d_1    + zeta4/8*(-Ck2 +2*Ck3 - Ck4)...
                      + zeta3/6*(3*Ck2 - 4*Ck3 + Ck4) - .5*zeta2*Ck2...
                      + (Ck2 +10*Ck3 + Ck4)/12;
        THET(kstr_vec(j)+2:K-1,j) = (CONT(kstr_vec(j)+1:K-2,j)+10*CONT(kstr_vec(j)+2:K-1,j)+CONT(kstr_vec(j)+3:K,j))/12;
        THET(K,j)       = (13*CONT(K,j)+15*CONT(K-1,j)-5*CONT(K-2,j)+CONT(K-3,j))/48;

    end
    
    CONT = zeros(K,m_0);
    for k = 1:m_0
        ThetTemp = fft([THET(1:K,k);zeros(K,1)]);
        for j = 1:m_0
            p = ifft(BETA(:,j,k).*ThetTemp);
            CONT(:,j) = CONT(:,j)+ p(1:K);
        end
    end
    
end

%%%%////////////////////////////////////////////////////////
%%%% Interpolate to find price at v0
%%%%////////////////////////////////////////////////////////
k_0 = 2;  %k_0 will satisfy v(k_0) <= v0 < v(k_0 +1)
while v0 > v(k_0) && k_0 < m_0
    k_0 = k_0+1;
end
k_0 = k_0 - 1;

%%% Cubic spline
% k_int = [(k_0-2) (k_0-1) k_0 (k_0+1) (k_0+2)];
% v_int = [v(k_int(1)) v(k_int(2)) v(k_int(3)) v(k_int(4)) v(k_int(5))];
% Vals_int = [CONT(nnot,(k_int(1))) CONT(nnot,(k_int(2))) CONT(nnot,(k_int(3))) CONT(nnot,(k_int(4))) CONT(nnot,(k_int(5)))];
% Vals_int = max(Vals_int, Gs(nnot));
% 
% price = spline(v_int,Vals_int,v0);

%%% Linear Interpolation
%%% NOTE: we have assumed that are not in Case II (see paper.. this case occurs if 0 < ln(W/S_0) < Delta), otherwise we need 2 interpolations
Vals_Interp = [max(CONT(nnot,(k_0)),Gs(nnot)) max(CONT(nnot,(k_0+1)),Gs(nnot))];
price = Vals_Interp(1) + (Vals_Interp(2)-Vals_Interp(1))*(v0 - v(k_0))/(v(k_0+1)-v(k_0));
end

