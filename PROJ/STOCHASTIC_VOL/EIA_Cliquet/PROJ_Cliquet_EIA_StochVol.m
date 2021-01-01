function price = PROJ_Cliquet_EIA_StochVol( numeric_param, M,r,q,T,psi_J,model, modparam, contract,contractParams)
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
%
% contract: 1 = sum of local caps
%           2 = sum of local caps & floors
%           3 = cliquet: local & global caps & floors
%           4 = cliquet: local floor & cap, global floor, NO GLOBAL CAP  
%           5 = MPP: ie monthly point-to-point or Monthly Cap Sum (Bernard, Li)
%
% contractParams:
% 	K  : Strike/Notional
% 	C  : Local Cap
% 	CG : Global cap
% 	F  : Local Floor
% 	FG : Global Floor
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = numeric_param.N;
alph = numeric_param.alph;
m_0 = numeric_param.m_0;
gridMethod = numeric_param.gridMethod;
gamma = numeric_param.gamma;
varGridMult = numeric_param.gridMultParam;

dx   = 2*alph/(N-1); a  = 1/dx;  %NOTE: this is just an initial guess, it is changed for contracts with both cap and floor
dt   = T/M;

xmin = (1-N/2)*dx; %Initialized value

%%% Contract Parameters (Not all of these apply to every contact type)
K  = contractParams.K;  % pricincal

C  = contractParams.C; %local cap
F  = contractParams.F; %local floor
CG = contractParams.CG; %global cap
FG = contractParams.FG; %global floor


lc = log(1 + C);
lf = log(1 + F);

%%% Choose xmin so that CAP lc is a member
klc = floor(a*(lc - xmin)) + 1;  %index of gridpoint left of CG
xklc = xmin + (klc - 1)*dx;
xmin = xmin + (lc - xklc); %Shift to the right (so index is the same)

klf = floor(a*(lf - xmin)) + 1;
%xklf = xmin + (klf - 1)*dx;  %NOTE: defined with the new xmin

if contract == 1 || contract == 5
    hlocalCF = @(x) (exp(x)- 1).*(x < lc) + C*(x>=lc);  %locally capped return
elseif contract == 2 || contract ==3 || contract ==4
    %NOTE: we should then possibly stretch the grid so that lf is a member
    if klc ~= klf
        dx = (lc - lf)/(klc - klf); a = 1/dx; 
        xmin = lf - (klf - 1)*dx;
    end   
    hlocalCF = @(x) F*(x<=lf) + (exp(x)- 1).*(x < lc).*(x > lf) + C*(x>=lc);  %locally capped and floored return   
end

A    = 32*a^4;
C_aN = A/N;

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
if contract == 2 || contract == 3 || contract == 4  
     leftGridPoint = lf-dx;  %this is the left bound of guassian quadrature grid 
     NNM = klc - klf + 3;  %this is N_Psi
    
elseif contract == 1 || contract == 5  %locally capped (but no floor)
    leftGridPoint = xmin; 
    NNM = klc+1;
else
    %NOTE: this can be made more efficient by putting an upper bound, to reflect lc
    leftGridPoint = xmin;
    NNM = N;  %% in this case, same dimension
end

    
PSI     = zeros(N-1,NNM);    %NOTE: we don't define the zeroth point in this code (used to be all ones)

%%%% Sample
Neta  = 5*(NNM) + 15;   %sample size
Neta5 = (NNM) + 3;
g2    = sqrt(5-2*sqrt(10/7))/6;
g3    = sqrt(5+2*sqrt(10/7))/6;
v1    = .5*128/225; v2 = .5*(322+13*sqrt(70))/900;  v3 = .5*(322 - 13*sqrt(70))/900;

thet                 = zeros(1,Neta);   %sample initialized
thet(5*(1:Neta5)-2)  = leftGridPoint -1.5*dx + dx*(0:Neta5-1);
thet(5*(1:Neta5)-4)  = leftGridPoint -1.5*dx + dx*(0:Neta5-1) - dx*g3;
thet(5*(1:Neta5)-3)  = leftGridPoint -1.5*dx + dx*(0:Neta5-1) - dx*g2;
thet(5*(1:Neta5)-1)  = leftGridPoint -1.5*dx + dx*(0:Neta5-1) + dx*g2;
thet(5*(1:Neta5))    = leftGridPoint -1.5*dx + dx*(0:Neta5-1) + dx*g3;

%%%% Weights
sig      = [-1.5-g3, -1.5-g2, -1.5, -1.5+g2, -1.5+g3, -.5-g3, -.5-g2, -.5, -.5+g2, -.5+g3,];
sig(1:5) = (sig(1:5) + 2).^3/6;
sig(6:10) = 2/3 - .5*(sig(6:10)).^3 - (sig(6:10)).^2;

sig([1 5 6 10]) = v3*sig([1 5 6 10]); sig([2 4 7 9]) = v2*sig([2 4 7 9]); sig([3 8]) = v1*sig([3 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%NEW STEP: multiple sig by Upsilon_{a,N}
sig = C_aN*sig;  %Note: I tested: the error of putting this here rather than directly after expansion is e-16 even in N = 2^12 case.. here is cheaper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fill Matrix
%%%% NOTE: this can be made MORE EFFICIENT by using symmetery of x^2

%zz  = exp(1i*dxi*log(1+exp(thet)));
%zz  = exp(1i*dxi*thet.^2); %% in general, 1i*dxh(thet)
zz = exp(1i*dxi*hlocalCF(thet));
thet   = zz; 

for j=1:N-1
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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Find phi_{Y_1}
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
xi    = dxi*(1:N-1)';  %% Redfined from above (doesnt have zeroth point)

b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
zeta  = (sin(xi/(2*a))./xi).^4./(b0 + b1*cos(xi/a) +b2*cos(2*xi/a) +b3*cos(3*xi/a));
hvec  = exp(-1i*xmin*xi).*zeta;  %dont define first element

PHIY_old = zeros(N-1,m_0);
PHIY_new = zeros(N-1,m_0);  %This serves intermediate storage for H,Beta,Phi_z
%BetaTemp = zeros(N,1);
PHI   = zeros(m_0,m_0,N-1);
grand = zeros(N-1,1); %%%temp vec

expFxi = exp(1i*F*xi);
expCxi = exp(1i*C*xi);

if contract == 2 || contract == 3 || contract == 4  
    for j = 1:m_0
        %Step 1: characteristic function of log return
        for n = 1:N-1
            PHIY_old(n,j) = sum(EXP_A(1:m_0,j,n+1));  %n+1 since EXP_A is defined for xi = 0
        end
        %Step 2: invert characteristic function of log return (ie this is beta)
        BetaTemp =  real(fft([1/A; PHIY_old(:,j).*hvec] )); %note: dont need PHIY_old(1,j)/A, since in this case PHIY_old(1,j) = 1

        %Step 3: Phi_{Y_1}^j 
        PHIY_new(:,j) = PSI*BetaTemp(klf-1:klc+1); 
        
        sumBetaLeft   = C_aN*sum(BetaTemp(1:klf-2));
        sumBetaRight  = 1 - sumBetaLeft - C_aN*sum(BetaTemp(klf-1:klc+1));
        PHIY_new(:,j) = PHIY_new(:,j) + expFxi*sumBetaLeft + expCxi*sumBetaRight;
    end
    
    % Define xiBig so that it can be added to a 3D matrix
    xiBigF = zeros(1,1,N-1);
    xiBigC = zeros(1,1,N-1);
    xiBigF(1,1,:) = expFxi;
    xiBigC(1,1,:) = expCxi;
    
    %%%NOTE: this next loop doesnt matter for M = 1
    if M > 1
        for j = 1:m_0
            for k=1:m_0
                %First Invert chf to get p_{j,k}
                for n=1:N-1  % MAKE this a .*
                    grand(n) = hvec(n)* EXP_A(k,j,n+1);  %n+1 since EXP_A has the first element defined
                end
                BetaTemp  = real(fft([EXP_A(k,j,1)/A; grand])); 

                PHI(j,k,:)   = PSI*BetaTemp(klf-1:klc+1); 
                sumBetaLeft  = C_aN*sum(BetaTemp(1:klf-2));
                sumBetaRight = C_aN*sum(BetaTemp(klc+2:N));
                %NOTE: can't do 1- here (see loop above, b/c the sum of beta is not 1 (its a conditional prob)
                PHI(j,k,:)   = PHI(j,k,:) + xiBigF*sumBetaLeft + xiBigC*sumBetaRight;        
            end
        end
    end
elseif contract == 5
    %%% ADD CODE
    fprintf('-------------------------------\n')
    fprintf('NOTE: HAVENT ADDED CODE FOR THIS CONTRACT\n\n\n')
    fprintf('-------------------------------\n')
end

%Main Recursion
for m =2:M
    for n = 1:N-1
        PHIY_new(n,:) = PHIY_new(n,:)*PHI(:,:,n).';
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Redfine ymin for the final inversion

%REDO FOR contract == 2 or ==3
if contract == 1 || contract ==2 %Sum of locally capped returns (ie sum min(C,Y_m)) 
    ymin = M*(exp((r-q)*dt) - 1) + (1-N/2)*dx;  %ie shifted by the expected sum of returns
    grid = ymin + dx*(0:N-1);
    
elseif contract == 3
    CminusF = CG - FG;
    ymin = FG - dx;
    kc = floor(a*(CG - ymin)) + 1;  %index of gridpoint left of CG
    z = a*(CG - (ymin +(kc -1)*dx));
    z2 = z^2; z3 = z*z2; z4 = z*z3; z5 = z*z4;
    
    theta = zeros(1,N/2);
    theta(1) = dx/120;
    theta(2) = dx*7/30;
    theta(3) = dx*121/120;
    theta(4:kc-2) = dx*(2:kc-4);
    k = kc-1;
    theta(k) = dx*(k*(-z4/24 + z3/6 - z2/4 +z/6 +23/24) - z5/30 +z4/6 - z3/3 +z2/3 - z/6 -59/30  )...
        + CminusF*(z-1)^4/24; 
    k = kc;
    theta(k) = dx*(k*(z4/8 - z3/3 +2*z/3 + .5)  +z5/10 -z4/2 +2*z3/3 +z2/3 -4*z/3 - 37/30)  ...
        + CminusF*(-z4/8 +z3/3 -2*z/3 +1/2);
    k = kc+1;
    theta(k) = dx*(k*(-z4/8 + z3/6 +z2/4 +z/6 +1/24) - z5/10 +z4/2 -z3/3 -2*z2/3 - z/2 - 2/15  )...
        + CminusF*(.5 +1/24*(3*z4 - 4*z3 - 6*z2 - 4*z + 11));
    k = kc+2;
    theta(k) = dx*(z5/30 +(k-4)*z4/24) + CminusF*(1-z4/24);
    theta(kc+3:N/2) = CminusF;
    
elseif contract == 4 || contract == 5  %NO GLOBAL CAP
    ymin = FG - dx; 
    theta = zeros(1,N/2);
    theta(1) = dx/120;
    theta(2) = dx*7/30;
    theta(3) = dx*121/120;
    theta(4:N/2) = dx*(2:N/2 - 2);
end


%%% Test with FILTER
applyFilter = 0;
if applyFilter == 1
    epsM = 1.2204e-16;   %matlabs machine epsilon
    alphaeps = -log(epsM);
    pp = 2; %order of the filter, must be even
    filter = exp(-alphaeps*(xi/(2*pi*a)).^pp);
    hvec = filter.*exp(-1i*ymin*xi).*zeta;
else
    hvec = exp(-1i*ymin*xi).*zeta;
end

% % %%%% Plot PHI %%%%%h = figure;
% h = figure;
% set(h,'defaultTextInterpreter','latex');
% dualGrand = .5*32*a^4*phi(2:N/1).*zeta(1:N/1-1);
% xGrid = [0; xi(1:N/1-1)];
% %plot(xGrid, abs([1; phi(2:N/1)]), xGrid, abs([1; dualGrand]), '--')
% p = plot(xGrid, abs([1; phi(2:N/1)]),'g', xGrid, abs([1; dualGrand]), 'r--', xGrid, abs([1; dualGrand.*filter]), 'b:');
% set(p,'MarkerSize',5,'LineWidth',1.1);
% ylabel('$|$Integrand$|$','Interpreter','LaTex')
% xlabel('$\xi$','Interpreter','LaTex')
% legend({'ChF','DualChF','Filtered'},'Interpreter','LaTex')
% %%%plot(abs(dualGrand.*filter))
% %%%plot(.5*32*a^4*zeta)


%BetaTemp = real(fft([1/A; hvec.*PHIY_new(:,1)])); 
% %%%%%%
% grid = ymin + dx*(0:N-1); %defined here only for plotting (grid is defined above in case of contract 1 or 2)
% plot(grid,C_aN*BetaTemp)
% %%%%%%


%%%%////////////////////////////////////////////////////////
%%%% Interpolate to find bracketing initial volatilities
%%%%////////////////////////////////////////////////////////
k_0 = 2;  %k_0 will satisfy v(k_0) <= v0 < v(k_0 +1)
while v0 >= v(k_0) && k_0 < m_0
    k_0 = k_0+1;
end
k_0 = k_0 - 1;
vals = [0 0];
ks = [k_0 k_0+1];

if contract == 1 || contract == 2
    BetaTemp = real(fft([1/A; hvec.*PHIY_new(:,ks(1))])); 
    vals(1) = K*exp(-r*T)*C_aN*grid(1:N)*BetaTemp(1:N);
    
    BetaTemp = real(fft([1/A; hvec.*PHIY_new(:,ks(2))])); 
    vals(2) = K*exp(-r*T)*C_aN*grid(1:N)*BetaTemp(1:N);
    
elseif contract == 3 || contract == 4 || contract == 5
    BetaTemp = real(fft([1/A; hvec.*PHIY_new(:,ks(1))])); 
    vals(1) = K*exp(-r*T)*(FG + C_aN*theta*BetaTemp(1:N/2));
    
    BetaTemp = real(fft([1/A; hvec.*PHIY_new(:,ks(2))])); 
    vals(2) = K*exp(-r*T)*(FG + C_aN*theta*BetaTemp(1:N/2));
end

%Linear Interpolation
if gridMethod == 5 || gridMethod == 6  %Then NO Interpolation Needed
    price = vals(1);
else
    price = vals(1) + (vals(2)-vals(1))*(v0 - v(k_0))/(v(k_0+1)-v(k_0));
end

end

