function price = MultiplicativeCliquet_PROJ(N, alph, M, r, T, rnCHF, contract, contractParams)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Cliquet-style options (Multiplicative Cliquets) using PROJ method
% Returns: price of contract
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Author: Justin Lars Kirkby
%
% ----------------------
% Contract/Model Params 
% ----------------------
% M = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
% r = interest rate (e.g. 0.05)
% T = time to maturity (in years, e.g. T=1)
% rnCHF = risk netural characteristic function (function handle with single argument)
% contract: 6 = Multiplicative Style Cliquet (e.g see Hieber)
% contractParams = container with the required params, such as cap and floor
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% N    = number of grid/basis points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
% alph = log-asset grid width param, grid with is 2*alph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx   = 2*alph/(N-1); a  = 1/dx;

xmin = (1-N/2)*dx; %Initialized value

%%% Contract Parameters (Not all of these apply to every contact type)
K  = contractParams.K;  % pricincal
C  = contractParams.C; %local cap
F  = contractParams.F; %local floor
Alpha = contractParams.Alpha;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

lc = log(C)/Alpha;
lf = log(F)/Alpha;

%%% Choose xmin so that CAP lc is a member
klc = floor(a*(lc - xmin)) + 1;  %index of gridpoint left of CG
xklc = xmin + (klc - 1)*dx;
xmin = xmin + (lc - xklc); %Shift to the right (so index is the same)

klf = floor(a*(lf - xmin)) + 1;

if contract == 6
    %NOTE: we should then possibly stretch the grid so that lf is a member
    if klc ~= klf
        dx = (lc - lf)/(klc - klf); a = 1/dx; 
        xmin = lf - (klf - 1)*dx;
    end   
    hlocalCF = @(x) F*(x<=lf) + exp(Alpha*x).*(x < lc).*(x > lf) + C*(x>=lc);  %locally capped and floored return   
end

A    = 32*a^4;
C_aN = A/N;
dxi  = 2*pi*a/N;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% PSI Matrix: 5-Point GAUSSIAN
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if contract == 6
    leftGridPoint = lf-dx;  %this is the left bound of guassian quadrature grid 
    NNM = klc - klf + 3;  %this is N_Psi
end

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
%%%% Fill Matrix
%%%% NOTE: this can be made MORE EFFICIENT by using symmetery of x^2

thet   = hlocalCF(thet); 


ThetaTilde =  sig(1)*(thet(1:5:Neta-19) + thet(20:5:Neta)) ...
          + sig(2)*(thet(2:5:Neta-18) + thet(19:5:Neta-1)) ...
          + sig(3)*(thet(3:5:Neta-17)  + thet(18:5:Neta-2)) ...
          + sig(4)*(thet(4:5:Neta-16)  + thet( 17:5:Neta-3)) ...
          + sig(5)*(thet(5:5:Neta-15)  + thet( 16:5:Neta-4)) ...
          + sig(6)*(thet(6:5:Neta-14)  + thet( 15:5:Neta-5)) ...
          + sig(7)*(thet(7:5:Neta-13)  + thet( 14:5:Neta-6)) ...
          + sig(8)*(thet(8:5:Neta-12)  + thet( 13:5:Neta-7)) ...
          + sig(9)*(thet(9:5:Neta-11)  + thet( 12:5:Neta-8)) ...
          + sig(10)*(thet(10:5:Neta-10)  + thet( 11:5:Neta-9));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi    = dxi*(1:N-1)';  %REDEFINED from above, excludes the zero at xi_1

b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
zeta  = (sin(xi/(2*a))./xi).^4./(b0 + b1*cos(xi/a) +b2*cos(2*xi/a) +b3*cos(3*xi/a));
hvec  = exp(-1i*xmin*xi).*zeta;  %dont define first element


AA = 1/A;
beta  = [AA; rnCHF(xi).*hvec];   %grand(end)=.5*grand(end);
beta  = real(fft(beta));


if contract == 6
    theta = zeros(1,N);
    theta(1:klf-2) = F;
    theta(klf-1:klc+1) = ThetaTilde;
    theta(klc+2:N) = C;
    
    price = theta*beta(1:N); %theta is 1xN/2
    price = K*exp(-r*T)*(C_aN*price)^M;
end

end

