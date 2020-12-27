function price = PROJ_DiscreteVariance_StochVol( numeric_param, M, r, T, K, psi_J, model, modparam, contract)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Discrete variance swaps and options using CTMC Approximation + PROJ method
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
% K  : strike  (used instead of K)
%
% contract: 
%           1 = Variance Swap, 
%           2 = Volatility Swap, 
%           3 = Call on Variance, 
%           4 = Put on Variance
%
% NOTE: right now only 1 and 3 are supported
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = numeric_param.N;
alph = numeric_param.alph;
m_0 = numeric_param.m_0;
gridMethod = numeric_param.gridMethod;
gamma = numeric_param.gamma;
varGridMult = numeric_param.gridMultParam;

dx   = 2*alph/(N-1);
dt   = T/M;
a    = 1/dx;
A    = 32*a^4;
C_aN = A/N;
xmin = (1-N/2)*dx; 

%%%%////////////////////////////////////////////////////////
%%%% Intialize Q matrix and variance set
%%%%////////////////////////////////////////////////////////
t = T/2;
[lx, v0, ux] = get_variance_grid_boundaries( model, modparam, t, gamma);

[ mu_func,  sig_func] = get_SV_variance_grid_diffusion_funcs( model,  modparam);
boundaryMethod = 1;
center = v0; %this is where grid clusters... we can experiment with other choices.. 

[Q,v]  = Q_Matrix_AllForms(m_0,mu_func,sig_func,lx,ux,gridMethod, varGridMult, center, boundaryMethod);


%%%%////////////////////////////////////////////////////////
%%%% Populate the Matrix Exponentials
%%%%////////////////////////////////////////////////////////
dxi    = 2*pi*a/N;
xi     = dxi*(0:N-1)';

[v1, v2, fv] = get_SV_matrix_expo_inputs( model,  modparam, psi_J, dt, v, dxi, r);
% Compute Matrix Exponentials for each xi(j)
EXP_A = get_SV_matrix_exponential( Q, dt, xi, v1, v2, fv, psi_J, m_0, N );


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% PSI Matrix: 5-Point GAUSSIAN
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
NNM = N;  %% in this case, same dimension
PSI = zeros(N-1,NNM);    %The first row (left undefined) will remain ones (since the chf will always be one at that point)

%%%% Sample
Neta  = 5*(NNM) + 15;   %sample size
Neta5 = (NNM) + 3;
g2    = sqrt(5-2*sqrt(10/7))/6;
g3    = sqrt(5+2*sqrt(10/7))/6;
v1    = .5*128/225; v2 = .5*(322+13*sqrt(70))/900;  v3 = .5*(322 - 13*sqrt(70))/900;

thet                 = zeros(1,Neta);   %sample initialized
thet(5*(1:Neta5)-2)  = xmin -1.5*dx + dx*(0:Neta5-1);
thet(5*(1:Neta5)-4)  = xmin -1.5*dx + dx*(0:Neta5-1) - dx*g3;
thet(5*(1:Neta5)-3)  = xmin -1.5*dx + dx*(0:Neta5-1) - dx*g2;
thet(5*(1:Neta5)-1)  = xmin -1.5*dx + dx*(0:Neta5-1) + dx*g2;
thet(5*(1:Neta5))    = xmin -1.5*dx + dx*(0:Neta5-1) + dx*g3;

%%%% Weights
sig      = [-1.5-g3, -1.5-g2, -1.5, -1.5+g2, -1.5+g3, -.5-g3, -.5-g2, -.5, -.5+g2, -.5+g3,];
sig(1:5) = (sig(1:5) + 2).^3/6;
sig(6:10) = 2/3 - .5*(sig(6:10)).^3 - (sig(6:10)).^2;

sig([1 5 6 10]) = v3*sig([1 5 6 10]); sig([2 4 7 9]) = v2*sig([2 4 7 9]); sig([3 8]) = v1*sig([3 8]);
sig = C_aN*sig;  

%%%% Fill Matrix  (NOTE: this can be made MORE EFFICIENT by using symmetery of x^2)
zz  = exp(1i*dxi*thet.^2); %% in general, 1i*dxh(thet)
thet   = zz; 

for j=1:N-1  %Note: first row is not ones anymore
    PSI(j,:) =  sig(1)*(thet(1:5:Neta-19) + thet(20:5:Neta)) ...
              + sig(2)*(thet(2:5:Neta-18) + thet(19:5:Neta-1)) ...
              + sig(3)*(thet(3:5:Neta-17)  + thet(18:5:Neta-2)) ...
              + sig(4)*(thet(4:5:Neta-16)  + thet( 17:5:Neta-3)) ...
              + sig(5)*(thet(5:5:Neta-15)  + thet( 16:5:Neta-4)) ...
              + sig(6)*(thet(6:5:Neta-14)  + thet( 15:5:Neta-5)) ...
              + sig(7)*(thet(7:5:Neta-13)  + thet( 14:5:Neta-6)) ...
              + sig(8)*(thet(8:5:Neta-12)  + thet( 13:5:Neta-7)) ...
              + sig(9)*(thet(9:5:Neta-11)  + thet( 12:5:Neta-8)) ...
              + sig(10)*(thet(10:5:Neta-10)  + thet( 11:5:Neta-9));

    thet = thet.*zz;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find phi_{Y_1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

xi    = dxi*(1:N-1)';  %REDEFINED from above, excludes the zero at xi_1
b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
zeta  = (sin(xi/(2*a))./xi).^4./(b0 + b1*cos(xi/a) +b2*cos(2*xi/a) +b3*cos(3*xi/a));
hvec  = exp(-1i*xmin*xi).*zeta;  %dont define first element

%Note: PHIY is not defined for xi = 0
PHIY_old = zeros(N-1,m_0);
PHIY_new = zeros(N-1,m_0);  %This serves intermediate storage for H,Beta,Phi_z

%Find beta in first stage (using chf of Phi_Y1)
for j = 1:m_0
   %Step 1: characteristic function of log return
   for n = 1:N-1
    PHIY_old(n,j) = sum(EXP_A(1:m_0,j,n+1));  %n+1 since EXP_A is defined for xi = 0
   end
   
   %Step 2: invert characteristic function of log return (ie this is beta)
   BetaTemp =  real(fft([1/A; PHIY_old(:,j).*hvec] )); %note: dont need PHIY_old(1,j)/A, since in this case PHIY_old(1,j) = 1
   
   %Step 3: Phi_{Y_1}^j 
   PHIY_new(:,j) = PSI*BetaTemp; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find PHI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
PHI  = ones(m_0,m_0,N-1); 
%%% use PHIY_old for temp storage
grand = zeros(N-1,1); %%%temp vec
for j = 1:m_0
    for k=1:m_0
        %First Invert chf to get p_{j,k}
        for n=1:N-1  % MAKE this a .*
            grand(n) = hvec(n)* EXP_A(k,j,n+1);  %n+1 since EXP_A has the first element defined
        end
        BetaTemp  = real(fft([EXP_A(k,j,1)/A; grand])); 
        PHI(j,k,:) = PSI*BetaTemp;
        %PHI(j,k,:) = PSI*real(fft([EXP_A(k,j,1)/A; grand]));
    end
end
clear EXP_A

for m =2:M
    for n = 1:N-1
        PHIY_new(n,:) = PHIY_new(n,:)*PHI(:,:,n).';
    end
end

%%%%////////////////////////////////////////////////////////
%%%% Interpolate to find bracketing initial volatilities
%%%%////////////////////////////////////////////////////////
k_0 = 2;  %k_0 will satisfy v(k_0) <= v0 < v(k_0 +1)
while v0 >= v(k_0) && k_0 < m_0
    k_0 = k_0+1;
end
k_0 = k_0 - 1;

%%%%%%%%%%%%%%%%%%%%%%%%% 
cubicTerminal = 1;  %Set to 1 to use cubic projection in final stage, else does linear 
%%%%%%%%%%%%%%%%%%%%%%%%% 
if contract == 1 %Variance Swap
    if cubicTerminal == 1
        xmin = -dx;
    else
        xmin = 0;
    end
elseif contract == 3 %Call on Variance
    if cubicTerminal == 1
        xmin = K*T - dx;
    else
        xmin = K*T;
    end
end

if cubicTerminal == 1
    if contract == 1 || contract == 3 %Variance Swap or Call on Variance   
        grid = -dx + dx*(0:N-1);  %NOTE: this holds for call on variance too, since we subtract out the KT

        grid(1) = grid(1)/24 + dx/20;
        grid(2) = dx*7/30;  %note: grid(2) = xmin + dx = 0;
        grid(3) = grid(3)*23/24 + dx/20;      
    end    
else  %Use Linear Projection at the end
    if contract == 1 || contract == 3 %Variance Swap or Call on Variance   
        grid = dx*(0:N-1);  %NOTE: this holds for call on variance too, since we subtract out the KT
        grid(1) = dx/6;
    end
    A = 24*a^2;
    C_aN  = A/N; 
    zeta = (sin(xi/(2*a))./xi).^2./(2+cos(xi/a));
end

if xmin ~= 0
    hvec = exp(-1i*xmin*xi).*zeta;
else
    hvec = zeta;
end

    
vals = [0 0];
ks = [k_0 k_0+1];    

for l = 1:2
    j = ks(l);
    BetaTemp  = real(fft([1/A; hvec.*PHIY_new(:,j)])); 
    if contract == 1 || contract == 3 %Variance Swap or call
        vals(l) = grid(1:N/2)*BetaTemp(1:N/2);
        vals(l) = C_aN*vals(l);
    end
end


if gridMethod == 5 || gridMethod == 6  %Then NO Interpolation Needed
    Approx = vals(1);
else
    Approx = vals(1) + (vals(2)-vals(1))*(v0 - v(k_0))/(v(k_0+1)-v(k_0));
end
    
if contract == 1 %Variance Swap
    price = Approx/T;  %can generalize to include anualization factor if needed
elseif contract == 3 %Variance Call
    price = exp(-r*T)*Approx/T;
else
    fprintf('Only contract types 1 and 3 are currently supported \n')
end

end
