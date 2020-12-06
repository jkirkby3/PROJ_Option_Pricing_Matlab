function Val = PROJ_Asian(N, alph, S_0, M, W, call, T, r, q, phiR, ER)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Arithmetic Asian Options using PROJ method
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% ----------------------
% Contract/Model Params 
% ----------------------
% S_0 = initial stock price (e.g. 100)
% W   = strike  (e.g. 100)
% r   = interest rate (e.g. 0.05)
% q   = dividend yield (e.g. 0.05)
% T   = time remaining until maturity (in years, e.g. T=1)
% M   = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
% call = 1 for call (else put)
% phiR =  risk neutral density of log return over time step dt = 1/M (function handle with single argument)
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% ER    = risk neutral expected return over increment dt=1/M (can set to zero if it is unknown)
% alph  = grid with is 2*alph
% N     = number of grid/basis points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = T/M;
dx = 2*alph/(N-1); a = 1/dx;

A    = 32*a^4;
C_aN = A/N;

%%% SHIFTS
x1    = zeros(1,M);
x1(1) = ER;
for m=2:M
    x1(m) = ER + log(1+exp(x1(m-1)));  %%BENHAMOU SHIFT
    %x1(m) = log(m) + .5*(m+1)*ER;    %% LOWER BOUND SHIFT derived in APROJ paper
end

Nm   = floor(a*(x1-ER));
x1   = ER + (1-N/2)*dx + Nm*dx;
NNM  = N + Nm(M-1);   %Number of columns of PSI

ystar = log((M+1)*W/S_0 -1);
nbar  = floor((ystar-x1(M))*a+1);

if nbar + 1 > N  
    % In this case, we fall off the grid when doing final integration, 
    % so you are deep ITM for the put, try to increase alph
    alph = 1.25*alph;
    Val = PROJ_Asian(N, alph, S_0, M, W, call, T, r, q, phiR, ER);
    return;
end

dxi   = 2*pi*a/N;
xi    = dxi*(1:(N-1))';
PhiR = [1; phiR(xi)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PSI Matrix: 5-Point GAUSSIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI     = zeros(N,NNM);    %The first row will remain ones
PSI(1,:) = ones(1,NNM);

%%%% Sample
Neta  = 5*(NNM) + 15;   %sample size
Neta5 = (NNM) + 3;
g2    = sqrt(5-2*sqrt(10/7))/6;
g3    = sqrt(5+2*sqrt(10/7))/6;
v1    = .5*128/225; 
v2    = .5*(322+13*sqrt(70))/900; 
v3    = .5*(322 - 13*sqrt(70))/900;


thet                 = zeros(1,Neta);   %sample initialized
thet(5*(1:Neta5)-2)  = x1(1) -1.5*dx + dx*(0:Neta5-1);
thet(5*(1:Neta5)-4)  = x1(1) -1.5*dx + dx*(0:Neta5-1) - dx*g3;
thet(5*(1:Neta5)-3)  = x1(1) -1.5*dx + dx*(0:Neta5-1) - dx*g2;
thet(5*(1:Neta5)-1)  = x1(1) -1.5*dx + dx*(0:Neta5-1) + dx*g2;
thet(5*(1:Neta5))    = x1(1) -1.5*dx + dx*(0:Neta5-1) + dx*g3;


%%%% Weights
sig      = [-1.5-g3, -1.5-g2, -1.5, -1.5+g2, -1.5+g3, -.5-g3, -.5-g2, -.5, -.5+g2, -.5+g3,];
sig(1:5) = (sig(1:5) + 2).^3/6;
sig(6:10) = 2/3 - .5*(sig(6:10)).^3 - (sig(6:10)).^2;

sig([1 5 6 10]) = v3*sig([1 5 6 10]); 
sig([2 4 7 9]) = v2*sig([2 4 7 9]); 
sig([3 8]) = v1*sig([3 8]);

%%%% Fill Matrix
zz  = exp(1i*dxi*log(1+exp(thet)));
thet   = zz; 


for j=2:N-1
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

%%-------------------------------------------------------------------------
b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
zeta  = (sin(xi/(2*a))./xi).^4./(b0 + b1*cos(xi/a) +b2*cos(2*xi/a) +b3*cos(3*xi/a));

AA = 1/A;
beta  = [AA; zeta.*PhiR(2:N).*exp(-1i*x1(1)*xi)];   %grand(end)=.5*grand(end);
beta  = real(fft(beta));

PhiR  = C_aN*PhiR;
beta  = PSI(:,1:N)*beta.*PhiR;  %Nm(1)=0

%%%%% Loop to find PSI_M
for m=3:M
    beta(2:N) = zeta.*beta(2:N).*exp(-1i*x1(m-1)*xi); beta(1) = AA;
    beta      = real(fft(beta));
    beta      = PSI(:,Nm(m-1)+1:Nm(m-1)+N)*beta.*PhiR;
end

%%-------------------------------------------------------------------------
%%%%% FINAL VALUE
C     = S_0/(M+1);
D     = W - C;
x1(M) = ystar- (nbar-1)*dx;

beta(2:N) = zeta.*beta(2:N).*exp(-1i*x1(M)*xi); beta(1)=AA;
beta      = real(fft(beta));

%%-------------------------------------------------------------------------
Cc1 = C*( exp(-7/4*dx)/54 + exp(-1.5*dx)/18 + exp(-1.25*dx)/2 + 7*exp(-dx)/27 )/20;

Cc2 = C*.05*(28/27 + exp(-7/4*dx)/54 + exp(-1.5*dx)/18 + exp(-1.25*dx)/2 + 14*exp(-dx)/27 ...
              + 121/54*exp(-.75*dx) + 23/18*exp(-.5*dx) + 235/54*exp(-.25*dx));
          
Cc3 = C*( (28 + 7*exp(-dx))/3 ...
              + ( 14*exp(dx) + exp(-7/4*dx) + 242*cosh(.75*dx) + 470*cosh(.25*dx))/12 ...
              +.25*(exp(-1.5*dx) + 9*exp(-1.25*dx) + 46*cosh(.5*dx)))/90;
          
Cc4 = C*( 14/3*(2+cosh(dx)) ...
              + .5*(cosh(1.5*dx) + 9*cosh(1.25*dx) +23*cosh(.5*dx))...
              +  1/6*(cosh(7/4*dx) + 121*cosh(.75*dx) +235*cosh(.25*dx)))/90;

%%-------------------------------------------------------------------------
G           = zeros(nbar+1,1);
E           = exp(ystar-(nbar-1)*dx+dx*(0:nbar));

G(nbar+1)   = D/24    - Cc1*E(nbar+1);
G(nbar)     = .5*D    - Cc2*E(nbar);
G(nbar-1)   = 23*D/24 - Cc3*E(nbar-1);
G(1:nbar-2) = D       - Cc4*E(1:nbar-2);
%%-------------------------------------------------------------------------

Val = C_aN*exp(-r*T)*sum(beta(1:nbar+1).*G);
if call==1  %Call Option, use Put-Call-Parity
    if r - q == 0
        mult = M + 1;
    else
        mult = (exp((r-q)*T*(1+1/M))-1)/(exp((r-q)*dt)-1);
    end
    Val = Val + C*exp(-r*T)*mult - W*exp(-r*T);
    
end
Val = max(0, Val);

end

