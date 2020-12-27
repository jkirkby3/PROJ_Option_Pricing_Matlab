function price = PROJ_Barrier_StochVol(numeric_param, call, down, S_0, W, H, M, r, T, psi_J, model, modparam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Barrier Options using CTMC Approximation + PROJ method
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
% H  : Barrier (either up and out, or down and out... No double barrier yet
% call : 1 for a call, else a put (easily can add digitals, etc)
% down : 1 for down-and-out, else up and out
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

N = numeric_param.N;
alph = numeric_param.alph;
m_0 = numeric_param.m_0;
gridMethod = numeric_param.gridMethod;
gamma = numeric_param.gamma;
gridMultParam = numeric_param.gridMultParam;

K    = N/2;
dx   = 2*alph/(N-1); a = 1/dx;
lws  = log(W/S_0);
dt   = T/M;

%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;

%%%% PAYOFF CONSTANTS-----------------------------------
varthet_01 = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
varthet_m10 = exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4)/18;
varthet_star = varthet_01 + varthet_m10;

%%%%////////////////////////////////////////////////////////
%%%% INITIALIZE LOG ASSET GRIDS (value grids)
%%%%////////////////////////////////////////////////////////
if down ==1  %%%%  DOWN & OUT 
    l    = log(H/S_0);  
    xmin = l;
    nnot = floor(1-xmin*a);   %index to left of x_0=0
    dx   = l/(1-nnot);  a    = 1/dx;
else %%% UP & OUT
    lws  = log(W/S_0);
    u    = log(H/S_0);
    nnot = floor(K-a*u);
    dx   = u/(K-nnot);
    a    = 1/dx;  
    xmin = u - (K-1)*dx;
end

nbar = floor(a*(lws - xmin)+1);
rho   = lws - (xmin+(nbar - 1)*dx);
zeta  = a*rho;

%%%%////////////////////////////////////////////////////////
%%%% Intialize Q matrix and variance set
%%%%////////////////////////////////////////////////////////

t = T/2;
[lx, v0, ux] = get_variance_grid_boundaries( model, modparam, t, gamma);

[ mu_func,  sig_func] = get_SV_variance_grid_diffusion_funcs( model,  modparam);
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

xi   = dxi*(1:N-1); %REDFINED FROM ABOVE
hvec = exp(-1i*zmin*xi).*(sin(xi/(2*a))./xi).^2./(2+cos(xi/a));

BETA  = zeros(N,m_0,m_0);  %to access the (j,k)th toeplitz data, use BETA(:,j,k)
grand = zeros(1,N-1);
      
%%% NOTE the (k,j) rather than (j,k)
for j=1:m_0
    for k = 1:m_0
        for n=1:N-1
            grand(n) = hvec(n)* EXP_A(k,j,n+1);  %n+1 since EXP_A has the first element defined
        end
        beta  = Cons2*real(fft([EXP_A(k,j,1)/(24*a2) grand]));   %%%%  NOTE: all toep matrices incorporate exp(-r*dt)
        toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];   
        BETA(:,j,k) = fft(toepM);
    end
end

%%%%////////////////////////////////////////////////////////
%%%% Initialiaze THETA (based on terminal payoff)
%%%%////////////////////////////////////////////////////////
if down ==1 %%% Down-and-out
    if call == 1  %%% Down-and-out CALL
        sigma = 1 - zeta; sigma_plus = (q_plus-.5)*sigma; sigma_minus = (q_minus-.5)*sigma;
        es1 = exp(dx*sigma_plus); es2 = exp(dx*sigma_minus);

        dbar_0 = .5 + zeta*(.5*zeta-1);
        dbar_1 = sigma*(1 - .5*sigma);
        d_0    = exp((rho+dx)*.5)*sigma^2/18*(5*((1-q_minus)*es2 +(1-q_plus)*es1) + 4);
        d_1    = exp((rho-dx)*.5)*sigma/18*(5*( (.5*(zeta+1) +sigma_minus)*es2 + (.5*(zeta+1) +sigma_plus)*es1 ) + 4*(zeta+1) );

        %%% NOTE: in first stage all columns are equal, so store them all in the THETA(1,:) position
        THET  = zeros(K,m_0);  %Each column corresponds to a regime state
        THET(nbar, 1)          = W*(exp(-rho)*d_0 - dbar_0);
        THET(nbar + 1, 1)      = W*(exp(dx-rho)*(varthet_01 +d_1) -(.5 + dbar_1) );
        THET(nbar + 2:K, 1)    = exp(xmin +dx*(nbar+1:K-1))*S_0*varthet_star - W;
    
    else %%% Down-and-out PUT   
        zeta_plus = zeta*q_plus; zeta_minus = zeta*q_minus;
        rho_plus = rho*q_plus; rho_minus = rho*q_minus;      
        ed1 = exp(rho_minus); ed2 = exp(rho/2); ed3 = exp(rho_plus);

        dbar_1 = zeta^2/2;
        dbar_0 = zeta - dbar_1;         %  dbar_1 = zeta + .5*((zeta - 1)^2 - 1);
        d_0    = zeta*(5*( (1-zeta_minus)*ed1 + (1-zeta_plus)*ed3 ) + 4*(2-zeta)*ed2)/18;
        d_1    = zeta*( 5*(zeta_minus*ed1 + zeta_plus*ed3) + 4*zeta*ed2 )/18;            

        %%% NOTE: in first stage all columns are equal, so store them all in the THETA(1,:) position
        THET  = zeros(K,m_0);  %Each column corresponds to a regime state
        THET(1, 1)        =  W/2 - H*varthet_01;
        THET(2:nbar-1, 1) =  W - exp(xmin +dx*(1:nbar-2))*S_0*varthet_star;
        THET(nbar, 1)     =  W*(.5 + dbar_0 - exp(-rho)*(varthet_m10 + d_0));
        THET(nbar + 1, 1) =  W*(dbar_1 - exp(- rho)*d_1);
    end
else % up-and-out
    if call == 1  %%% Up-and-out CALL 
        sigma = 1 - zeta; sigma_plus = (q_plus-.5)*sigma; sigma_minus = (q_minus-.5)*sigma;
        es1 = exp(dx*sigma_plus); es2 = exp(dx*sigma_minus);

        dbar_0 = .5 + zeta*(.5*zeta-1);
        dbar_1 = sigma*(1 - .5*sigma);
        d_0 = exp((rho+dx)*.5)*sigma^2/18*(5*((1-q_minus)*es2 +(1-q_plus)*es1) + 4);
        d_1 = exp((rho-dx)*.5)*sigma/18*(5*( (.5*(zeta+1) +sigma_minus)*es2 + (.5*(zeta+1) +sigma_plus)*es1 ) + 4*(zeta+1) );

        %%% NOTE: in first stage all columns are equal, so store them all in the THETA(1,:) position
        THET  = zeros(K,m_0);  %Each column corresponds to a regime state
        THET(nbar, 1)          = W*(exp(-rho)*d_0 - dbar_0);
        THET(nbar + 1, 1)      = W*(exp(dx-rho)*(varthet_01 +d_1) -(.5 + dbar_1) );
        THET(nbar + 2:K-1, 1)  = exp(xmin +dx*(nbar+1:K-2))*S_0*varthet_star - W;
        THET(K, 1)             = H*varthet_m10 - .5*W;

    else %%% Up-and-out PUT 
        zeta_plus = zeta*q_plus; zeta_minus = zeta*q_minus;
        rho_plus = rho*q_plus; rho_minus = rho*q_minus;
        ed1 = exp(rho_minus); ed2 = exp(rho/2); ed3 = exp(rho_plus);

        dbar_1 = zeta^2/2;
        dbar_0 = zeta - dbar_1;         %  dbar_1 = zeta + .5*((zeta - 1)^2 - 1);
        d_0    = zeta*(5*( (1-zeta_minus)*ed1 + (1-zeta_plus)*ed3 ) + 4*(2-zeta)*ed2)/18;
        d_1    = zeta*( 5*(zeta_minus*ed1 + zeta_plus*ed3) + 4*zeta*ed2 )/18;            

        %%% NOTE: in first stage all columns are equal, so store them all in the THETA(1,:) position
        THET  = zeros(K,m_0);  %Each column corresponds to a regime state
        THET(1:nbar-1, 1) =  W - exp(xmin +dx*(0:nbar-2))*S_0*varthet_star;
        THET(nbar, 1)     =  W*(.5 + dbar_0 - exp(-rho)*(varthet_m10 + d_0));
        THET(nbar + 1, 1) =  W*(dbar_1 - exp(- rho)*d_1);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize Continuation Value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CONT = zeros(K,m_0);
%%% Note: First column of THET used to store ThetM at intialization
ThetTemp = fft([THET(1:K, 1);zeros(K,1)]);
for j=1:m_0
    for k = 1:m_0
        p = ifft(BETA(:,j,k).*ThetTemp); 
        CONT(:,j) = CONT(:,j)+ p(1:K);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOOP through time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = M-2:-1:0
    %Step 1: update THETA
    for j=1:m_0
        THET(1,j) = (13*CONT(1,j) + 15*CONT(2,j) - 5*CONT(3,j) + CONT(4,j))/48;
        THET(K,j) = (13*CONT(K,j) + 15*CONT(K-1,j) - 5*CONT(K-2,j) + CONT(K-3,j))/48;
        THET(2:K-1,j) = (CONT(1:K-2,j) +10*CONT(2:K-1,j) + CONT(3:K,j))/12;
    end
    
    %Step 2: sum up the convolutions
    CONT = zeros(K,m_0);
    for k = 1:m_0
        ThetTemp = fft([THET(1:K,k);zeros(K,1)]);
        for j = 1:m_0
            p = ifft(BETA(:,j,k).*ThetTemp);
            CONT(:,j) = CONT(:,j)+ p(1:K);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate to find price at v0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_0 = 2;  %k_0 will satisfy v(k_0) <= v0 < v(k_0 +1)
while v0 > v(k_0) && k_0 < m_0
    k_0 = k_0+1;
end
k_0 = k_0 - 1;

%%% Cubic Interpolation
% k_int = [(k_0-1) k_0 (k_0+1)];
% v_int = [v(k_int(1)) v(k_int(2)) v(k_int(3))];
% Vals_int = [CONT(nnot,(k_int(1))) CONT(nnot,(k_int(2))) CONT(nnot,(k_int(3)))];
% Val = spline(v_int,Vals_int,v0);

%%% Linear Interpolation (w.r.t initial regime)
Vals_Interp = [CONT(nnot,(k_0)) CONT(nnot,(k_0+1))];  %note: nnot is grid point except in special case, which we assume doesn't occur (for now)
price = Vals_Interp(1) + (Vals_Interp(2)-Vals_Interp(1))*(v0 - v(k_0))/(v(k_0+1)-v(k_0));


end



  