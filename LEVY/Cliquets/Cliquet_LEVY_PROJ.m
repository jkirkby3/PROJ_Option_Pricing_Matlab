function price = Cliquet_LEVY_PROJ(N, alph, M, r, q, T, rnCHF, contract, contractParams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Cliquet-style options (Additive Cliquets) using PROJ method
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% ----------------------
% Contract/Model Params 
% ----------------------
% M = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
% r = interest rate (e.g. 0.05)
% q = dividend yield (e.g. 0.05)
% T = time to maturity (in years, e.g. T=1)
% rnCHF = risk netural characteristic function (function handle with single argument)
% contract: 1 = sum of local caps
%           2 = sum of local caps & floors
%           3 = cliquet: local & global caps & floors
%           4 = cliquet: local floor & cap, global floor, NO GLOBAL CAP (e.g. like Wilmott)  
%           5 = MPP: ie monthly point-to-point  or Monthly Cap Sum (Bernard, Li)
% contractParams - container with the required params, such as cap and floor
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% N = number of grid/basis points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
% alph = log-asset grid width param, grid with is 2*alph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx   = 2*alph/(N-1); a  = 1/dx;
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
dxi  = 2*pi*a/N;

    
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
sig = C_aN*sig;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fill Matrix
%%%% NOTE: this can be made MORE EFFICIENT by using symmetery of x^2

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

xi    = dxi*(1:N-1)';  %REDEFINED from above, excludes the zero at xi_1

b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
zeta  = (sin(xi/(2*a))./xi).^4./(b0 + b1*cos(xi/a) +b2*cos(2*xi/a) +b3*cos(3*xi/a));
hvec  = exp(-1i*xmin*xi).*zeta;  %dont define first element


AA = 1/A;
beta  = [AA; rnCHF(xi).*hvec];   %grand(end)=.5*grand(end);
beta  = real(fft(beta));

if contract == 2 || contract == 3 || contract == 4  
    phi = PSI*beta(klf-1:klc+1); 
    sumBetaLeft = C_aN*sum(beta(1:klf-2));
    sumBetaRight = 1 - sumBetaLeft - C_aN*sum(beta(klf-1:klc+1));
    phi = phi + exp(1i*F*xi)*sumBetaLeft + exp(1i*C*xi)*sumBetaRight;
elseif contract == 1 || contract == 5
    phi = PSI*beta(1:klc+1); 
    sumBetaRight = C_aN*sum(beta(klc+2:N));
    phi = phi + exp(1i*C*xi)*sumBetaRight;
else
    phi = PSI*beta; 
end


phi = phi.^M;  %NOTE: this is the entire convolution step, b/c indepdendent stationary increments


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Redfine xmin for the final inversion

%REDO FOR contract == 2 or ==3
if contract == 1 || contract ==2 %Sum of locally capped returns (ie sum min(C,Y_m)) 
    ymin = M*(exp((r-q)*dt) - 1) + (1-N/2)*dx;  %ie shifted by the expected sum of returns
    %grid = ymin + dx*(0:N-1);
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


%%% Filtering (optional)
applyFilter = 0;  %HARDCODED - can switch on/off for testing
if applyFilter == 1
    epsM = 1.2204e-16;   %matlabs machine epsilon
    alphaeps = -log(epsM);
    pp = 4; %order of the filter, must be even
    filter = exp(-alphaeps*(xi/(2*pi*a)).^pp);
    hvec = filter.*exp(-1i*ymin*xi).*zeta;
else
    hvec = exp(-1i*ymin*xi).*zeta;
end


beta = real(fft([1/A; hvec.*phi])); 

% Final pricing step, depends on contract type
if contract == 1 || contract == 2
    grid = ymin + dx*(0:N-1);
    price = grid(1:N)*beta(1:N);
    price = K*exp(-r*T)*C_aN*price;
elseif contract == 3 || contract == 4 || contract == 5
    price = theta*beta(1:N/2); %theta is 1xN/2
    price = K*exp(-r*T)*(FG + C_aN*price);
end


end

