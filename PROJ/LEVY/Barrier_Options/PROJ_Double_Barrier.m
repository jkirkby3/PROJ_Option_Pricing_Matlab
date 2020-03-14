function [price] = PROJ_Double_Barrier(N,alph,call,L,U, S_0,W,M,T,r,rnCHF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Discrete Double Barrier Options using PROJ method
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
% [L,U] = barriers
% rnCHF = risk netural characteristic function (function handle with single argument)
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% N     = number of grid/basis points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


l = log(L/S_0);  u = log(U/S_0);  dt = T/M;

K = N/2;
dx = (u-l)/(K-1);
a  = 1/dx;

E  = ceil(2*alph/(u-l));
if M<= 12
    E = min(E,4);
else
    E = min(E,3);
end


a2    = a^2;
N_Ee  = E*N;
Cons2 = 24*a2*exp(-r*dt)/N_Ee ;


grand = @(w)rnCHF(w).*(sin(w/(2*a))./w).^2./(2+cos(w/a));
dw    = 2*pi*a/N_Ee ;
omega = (dw: dw: (N_Ee -1)*dw);      %We calcuate coefficient of w=0 explicitly
zmin  = (1 - E*K)*dx;           %K corresponds to zero
beta  = Cons2*real(fft([1/(24*a^2) exp(-1i*zmin*omega).*feval(grand,omega)]));

toepM = [beta(E*K:-1:(E-1)*K+1)';0; beta((E+1)*K-1:-1:E*K+1)'];
toepM = fft(toepM);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Payoff Coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin = l;
nnot = floor(1-xmin*a);   %index to left of x_0=0
lws = log(W/S_0);
nbar = floor(a*(lws - xmin)+1);
rho   = lws - (xmin+(nbar - 1)*dx);
zeta  = a*rho;
xnbar = xmin + (nbar - 1)*dx;

Thet  = zeros(K,1);
Cons3 = 1/48;
Cons4 = 1/12;

%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;


%%%% PAYOFF CONSTANTS-----------------------------------
varthet_01 = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
varthet_m10 = exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4)/18;
varthet_star = varthet_01 + varthet_m10;

if call==1  %DBC
    sigma = 1 - zeta; sigma_plus = (q_plus-.5)*sigma; sigma_minus = (q_minus-.5)*sigma;
    es1 = exp(dx*sigma_plus); es2 = exp(dx*sigma_minus);
    dbar_0 = .5 + zeta*(.5*zeta-1);
    dbar_1 = sigma*(1 - .5*sigma);

    d_0 = exp((rho+dx)*.5)*sigma^2/18*(5*((1-q_minus)*es2 +(1-q_plus)*es1) + 4);
    d_1 = exp((rho-dx)*.5)*sigma/18*(5*( (.5*(zeta+1) +sigma_minus)*es2 + (.5*(zeta+1) +sigma_plus)*es1 ) + 4*(zeta+1) );
    
   
    Thet(nbar)         = W*(exp(-rho)*d_0 - dbar_0);
    Thet(nbar + 1)     = W*(exp(dx-rho)*(varthet_01 +d_1) -(.5 + dbar_1) );
    Thet(nbar + 2:K-1) = exp(xmin +dx*(nbar+1:K-2))*S_0*varthet_star - W;
    Thet(K)            = U*varthet_m10 - W/2;
    %%%%%%%
    p = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
    Val = p(1:K);
    %%%%%%%
    for m=M-2:-1:0
        Thet(1)      = Cons3*(13*Val(1)+15*Val(2)-5*Val(3)+Val(4));
        Thet(K)      = Cons3*(13*Val(K)+15*Val(K-1)-5*Val(K-2)+Val(K-3));
        Thet(2:K -1) = Cons4*(Val(1:K-2)+10*Val(2:K-1)+Val(3:K));
        %%%%%%%
        p        = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
        Val = p(1:K);
    end
    
    
else  %DBP
    zeta_plus = zeta*q_plus; zeta_minus = zeta*q_minus;
    rho_plus = rho*q_plus; rho_minus = rho*q_minus;

    ed1 = exp(rho_minus); ed2 = exp(rho/2); ed3 = exp(rho_plus);

    dbar_1 = zeta^2/2;
    dbar_0 = zeta - dbar_1;         %  dbar_1 = zeta + .5*((zeta - 1)^2 - 1);
    d_0    = zeta*(5*( (1-zeta_minus)*ed1 + (1-zeta_plus)*ed3 ) + 4*(2-zeta)*ed2)/18;
    d_1    = zeta*( 5*(zeta_minus*ed1 + zeta_plus*ed3) + 4*zeta*ed2 )/18;            

    Thet(1)        =  W/2 - L*varthet_01;
    Thet(2:nbar-1) =  W - exp(xmin +dx*(1:nbar-2))*S_0*varthet_star;
    Thet(nbar)     =  W*(.5 + dbar_0 - exp(-rho)*(varthet_m10 + d_0));
    Thet(nbar + 1) =  W*(dbar_1 - exp(- rho)*d_1);
    %%%%%%%
    p   = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
    Val = p(1:K);
    %%%%%%%
    for m=M-2:-1:0
        Thet(1)      = Cons3*(13*Val(1)+15*Val(2)-5*Val(3)+Val(4));
        Thet(K)      = Cons3*(13*Val(K)+15*Val(K-1)-5*Val(K-2)+Val(K-3));
        Thet(2:K -1) = Cons4*(Val(1:K-2)+10*Val(2:K-1)+Val(3:K));
        %%%%%%%
        p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
        Val  = p(1:K);
    end   
end

xnot = l+(nnot-1)*dx;
xs = [xnot-2*dx,xnot-dx,xnot,xnot+dx,xnot+2*dx];
ys = [Val(nnot-2),Val(nnot-1),Val(nnot),Val(nnot+1),Val(nnot+2)];
price = spline(xs,ys,0);


end

