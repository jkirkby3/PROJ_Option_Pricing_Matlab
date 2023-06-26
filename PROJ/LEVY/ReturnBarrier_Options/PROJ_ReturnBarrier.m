function [price, QM] = PROJ_ReturnBarrier(N, alpha_T, alpha_dt, M, r, q, T, S_0, rnCHF, contractParams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Return Barrier Options using PROJ method
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% Reference: (1) "The Return Barrier and Return Timer Option with Pricing Under
%              Levy Processes", J.L. Kirkby and J-P Aguilar, Expert Systems with Applications, 2023
%            (2) "Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform", J.L. Kirkby, SIAM J. Financial Math., 2015
% 
% ----------------------
% Contract/Model Params 
% ----------------------
% M = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
% r = interest rate (e.g. 0.05)
% q = dividend yield (e.g. 0.05)
% T = time to maturity (in years, e.g. T=1)
% rnCHF = risk netural characteristic function (function handle with single argument)
% contractParams - container with the required params
%                - l, u
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% N = number of grid/basis points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
% alph = log-asset grid width param, grid with is 2*alph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx   = 2*alpha_dt/(N-1); a  = 1/dx;
xmin = (1-N/2)*dx; %Initialized value

%%% Contract Parameters (Not all of these apply to every contact type)
K  = contractParams.K;  % strike of option
lf  = contractParams.F; %local floor
lc  = contractParams.C; %local cap

lf = max(xmin+2*dx, lf);
lc = min(lc, xmin + (N-3)*dx);

%%% Choose xmin so that lc and lf are members
klc = floor(a*(lc - xmin)) + 1;  %index of gridpoint left of C
xklc = xmin + (klc - 1)*dx;
xmin = xmin + (lc - xklc); %Shift to the right (so index is the same)

klf = floor(a*(lf - xmin)) + 1;

if klc ~= klf
    dx = (lc - lf)/(klc - klf); a = 1/dx; 
    xmin = lf - (klf - 1)*dx;
end  

A    = 32*a^4;
C_aN = A/N;
dxi  = 2*pi*a/N;

%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% PSI Matrix: 5-Point GAUSSIAN
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
leftGridPoint = lf-dx;  %this is the left bound of guassian quadrature grid 
NNM = klc - klf + 3;  %this is N_Psi
PSI     = zeros(N,NNM);

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
thet0 = thet;

%%%% Weights
sig      = [-1.5-g3, -1.5-g2, -1.5, -1.5+g2, -1.5+g3, -.5-g3, -.5-g2, -.5, -.5+g2, -.5+g3,];
sig(1:5) = (sig(1:5) + 2).^3/6;
sig(6:10) = 2/3 - .5*(sig(6:10)).^3 - (sig(6:10)).^2;

sig([1 5 6 10]) = v3*sig([1 5 6 10]); sig([2 4 7 9]) = v2*sig([2 4 7 9]); sig([3 8]) = v1*sig([3 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% multiple sig by Upsilon_{a,N}
sig = C_aN*sig;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fill Matrix
dx_T   = 2*alpha_T/(N-1); 
a_T  = 1/dx_T;
dxi_T  = 2*pi*a_T/N;

zz = exp(1i*dxi_T*thet);

thet   = ones(1,length(thet)); 
thet(thet0 <= lf) = 0;
thet(thet0 >= lc) = 0;

for j=1:N
    PSI(j,:) = FillQuadratureRow(sig, thet, Neta);
    thet = thet.*zz;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Find phi_{Y_1}
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
xi    = dxi*(1:N-1)';  
zeta = get_zeta(xi, a);

applyFilter = 1;  %HARDCODED - can switch on/off for testing
if applyFilter == 1
    epsM = 1.2204e-16;   %matlabs machine epsilon
    alphaeps = -log(epsM);
    pp = 6; %order of the filter, must be even
    filter = exp(-alphaeps*(xi/(2*pi*a)).^pp);
    hvec = filter.*exp(-1i*xmin*xi).*zeta;
else
    hvec  = exp(-1i*xmin*xi).*zeta;  %dont define first element
end

AA = 1/A;
beta  = [AA; rnCHF(xi).*hvec];   %grand(end)=.5*grand(end);
beta  = real(fft(beta));
  
phi = PSI*beta(klf-1:klc+1);  % Only one basis element to the left/right of each boundary overlaps the boundary
Q = real(phi(1));

phi = (phi).^M;  %NOTE: this is the entire convolution step, b/c indepdendent stationary increments
% QM = real(phi(1));

phirmi = 0;
if contractParams.call == 1
    % Compute Expected value of R|[l,u]
    thet =exp(thet0) ;
    thet(thet0 <= lf) = 0;
    thet(thet0 >= lc) = 0;
    PSI_mi =  FillQuadratureRow(sig, thet, Neta);
    phirmi = PSI_mi*beta(klf-1:klc+1);  % Only one basis element to the left/right of each boundary overlaps the boundary
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now apply cubic PROJ method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
price = final_price(dx_T, K, N, S_0, phi, Q, phirmi, M, r, T, contractParams.call);
end

function price = final_price(dx_T, K, N, S_0, phi, Q, phirmi, M, r, T,call)
a_T  = 1/dx_T;
dxi_T  = 2*pi*a_T/N;
A_T    = 32*a_T^4;

QM = Q.^M;

%%%%%%%%%%%%%%%%

dx = dx_T;
a = a_T;

% TODO: estimate c1
c1 = 0;
W = K;
lws = log(W/S_0);
lam = c1 -(N/2 -1)*dx;
nbar = floor(a*(lws-lam)+1);
if nbar>=N
    nbar = N-1;
end

ymin = lws - (nbar-1)*dx;
xi_T    = dxi_T*(1:N-1)';  %REDEFINED from above, excludes the zero at xi_1
zeta = get_zeta(xi_T, a);
hvec = exp(-1i*ymin*xi_T).*zeta;

beta = real(fft([QM/A_T; hvec.*phi(2:end)])); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------   CUBIC  ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

G = zeros(1,nbar +1); 
G(nbar +1) = W*(1/24 - 1/20*exp(dx)*(exp(-7/4*dx)/54 + exp(-1.5*dx)/18 + exp(-1.25*dx)/2 + 7*exp(-dx)/27));

G(nbar )   =  W*(.5 -.05*(28/27 + exp(-7/4*dx)/54 + exp(-1.5*dx)/18 + exp(-1.25*dx)/2 + 14*exp(-dx)/27 ...
            + 121/54*exp(-.75*dx) + 23/18*exp(-.5*dx) + 235/54*exp(-.25*dx)));

G(nbar -1) = W*( 23/24 - exp(-dx)/90*( (28 + 7*exp(-dx))/3 ...
            + ( 14*exp(dx) + exp(-7/4*dx) + 242*cosh(.75*dx) + 470*cosh(.25*dx))/12 ...
            +.25*(exp(-1.5*dx) + 9*exp(-1.25*dx) + 46*cosh(.5*dx))) );

G(1: nbar -2) = W - S_0*exp(ymin +dx*(0:nbar-3))/90*( 14/3*(2+cosh(dx)) ...
                + .5*(cosh(1.5*dx) + 9*cosh(1.25*dx) +23*cosh(.5*dx))...
                +  1/6*(cosh(7/4*dx) + 121*cosh(.75*dx) +235*cosh(.25*dx)));

Cons = 32*a^4;

disc = exp(-r*T);
price_put = Cons*disc/N*G*beta(1:length(G));

if call == 1  % Use put-call parity
    price = price_put + disc*S_0*phirmi.^M - disc*W*QM;
else
    price = price_put;
end

price = max(price, 0);  % Protect against deep out of money case
end

function zeta = get_zeta(xi, a)
b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
zeta  = (sin(xi/(2*a))./xi).^4./(b0 + b1*cos(xi/a) +b2*cos(2*xi/a) +b3*cos(3*xi/a));
end

function quadrature_row = FillQuadratureRow(sig, thet, Neta)
quadrature_row =  sig(1)*(thet(1:5:Neta-19) + thet(20:5:Neta)) ...
          + sig(2)*(thet(2:5:Neta-18) + thet(19:5:Neta-1)) ...
          + sig(3)*(thet(3:5:Neta-17)  + thet(18:5:Neta-2)) ...
          + sig(4)*(thet(4:5:Neta-16)  + thet( 17:5:Neta-3)) ...
          + sig(5)*(thet(5:5:Neta-15)  + thet( 16:5:Neta-4)) ...
          + sig(6)*(thet(6:5:Neta-14)  + thet( 15:5:Neta-5)) ...
          + sig(7)*(thet(7:5:Neta-13)  + thet( 14:5:Neta-6)) ...
          + sig(8)*(thet(8:5:Neta-12)  + thet( 13:5:Neta-7)) ...
          + sig(9)*(thet(9:5:Neta-11)  + thet( 12:5:Neta-8)) ...
          + sig(10)*(thet(10:5:Neta-10)  + thet( 11:5:Neta-9));
end


