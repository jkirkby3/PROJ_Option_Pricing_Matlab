function [prob,spread] = PROJ_CDS(R, L, M, T, r, N, alph, mult, rnCHF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Calc Fair Spread of Credit Default Swaps (and default probabilities) using PROJ method
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: prob = probability of default on [0,T]
%          spread = fair (par) CDS spread in basis points (ie *10000)
% Author: Justin Lars Kirkby
%
% ----------------------
% Contract/Model Params 
% ----------------------
% r     = interest rate (e.g. 0.05)
% T     = time remaining until maturity (in years, e.g. T=1)
% M     = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including initial time)
% R     = recovery rate, 0<R<1
% L     = percentage of initial firm value that leads to default, 0<L<1
% rnCHF = risk netural characteristic function (function handle with single argument)
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% alph  = grid with is 2*alph
% N     = number of grid/basis points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt   = T/M;
l    = log(L);

K = N/2;
dx = 2*alph/(N-1); a = 1/dx;

xmin = l;
nnot = floor(1-xmin*a);   %index to left of x_0=0
dx   = l/(1-nnot);
a    = 1/dx;

zmin  = (1 - K)*dx; 
Nmult = mult*N;
dw    = 2*pi*a/Nmult;
grand = dw*(1:Nmult-1);

Thet  = ones(K,1);  Thet(1) = .5;
b1    = 1/48;  b2    = 1/12; 

a2    = a^2;
Cons = 24*a2/Nmult;

grand = exp(-1i*zmin*grand).*rnCHF(grand).*(sin(grand/(2*a))./grand).^2./(2+cos(grand/a));
beta  = Cons*real(fft([1/(24*a2) grand]));  
toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];   toepM = fft(toepM);

Thetbar1 = cumsum(beta(2*K:-1:K +1))';

p = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
Val = p(1:K) + Thetbar1;

Probs    = zeros(1,M);
Probs(1) =  Val(nnot);

for m=M-2:-1:0
    Thet(1)      = b1*(13*Val(1)+15*Val(2)-5*Val(3)+Val(4));
    Thet(K)      = b1*(13*Val(K)+15*Val(K-1)-5*Val(K-2)+Val(K-3));
    Thet(2:K -1) = b2*(Val(1:K-2)+10*Val(2:K-1)+Val(3:K));
    %%%%%%%
    p        = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
    Val(1:K) = p(1:K)+ Thetbar1;
    Probs(M-m) = Val(nnot);
end

prob = Val(nnot);
denom = dt*(.5 + sum( exp(-r*dt*(1:M-1)).*Probs(1:M-1))  + .5*Probs(M)*exp(-r*T));  % trapezoid rule
spread = 10000*(1-R)*( (1 - exp(-r*T)*Probs(M))/denom - r);

end

