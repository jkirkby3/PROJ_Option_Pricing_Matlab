function price = PROJ_European(order, N, alph, r, q, T, S_0, W, call, rnCHF, c1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European Options using PROJ method
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% References:  1) Efficient Option Pricing by Frame Duality with the Fast Fourier Transform. 
%                 SIAM J. Financial Math (2015), Kirkby, J.L
%              2) Robust option pricing with characteristic functions and the B-spline order of density projection,
%                 J. Compuational Finance (2017), Kirkby, J.L
%
% ----------------------
% Contract/Model Params 
% ----------------------
% S_0 = initial stock price (e.g. 100)
% W   = strike  (e.g. 100)
% r   = interest rate (e.g. 0.05)
% q   = dividend yield (e.g. 0.05)
% T   = time remaining until maturity (in years, e.g. T=1)
% call  = 1 for call (else put)
% rnCHF = risk netural characteristic function (function handle with single argument)
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% order = 0,1,2,3  (Order of spline: Haar,Linear,Quadratic,Cubic
% c1    = mean return, first cumulant (incorporates T.. can safefly set to zero for small maturities, say T<=2)
% alph  = grid with is 2*alph
% N     = number of grid/basis points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = 2*alph/(N-1); a = 1/dx;

lws = log(W/S_0);
lam = c1 -(N/2 -1)*dx;
nbar = floor(a*(lws-lam)+1);
if nbar>=N
    nbar = N-1;
end
xmin = lws - (nbar-1)*dx;

dw = 2*pi/(N*dx);
omega = (dw: dw: (N-1)*dw);  %We calcuate coefficient of w=0 explicitly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------   CUBIC  ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if order ==3  
    b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
    grand = @(w)rnCHF(w).*(sin(w/(2*a))./w).^4./(b0 + b1*cos(w/a) +b2*cos(2*w/a) +b3*cos(3*w/a));
    beta  = real(fft([1/(32*a^4) exp(-1i*xmin*omega).*feval(grand,omega)]));
    
    G = zeros(1,nbar +1); 
    G(nbar +1) = W*(1/24 - 1/20*exp(dx)*(exp(-7/4*dx)/54 + exp(-1.5*dx)/18 + exp(-1.25*dx)/2 + 7*exp(-dx)/27));
    
    G(nbar )   =  W*(.5 -.05*(28/27 + exp(-7/4*dx)/54 + exp(-1.5*dx)/18 + exp(-1.25*dx)/2 + 14*exp(-dx)/27 ...
                + 121/54*exp(-.75*dx) + 23/18*exp(-.5*dx) + 235/54*exp(-.25*dx)));
            
    G(nbar -1) = W*( 23/24 - exp(-dx)/90*( (28 + 7*exp(-dx))/3 ...
                + ( 14*exp(dx) + exp(-7/4*dx) + 242*cosh(.75*dx) + 470*cosh(.25*dx))/12 ...
                +.25*(exp(-1.5*dx) + 9*exp(-1.25*dx) + 46*cosh(.5*dx))) );
            
    G(1: nbar -2) = W - S_0*exp(xmin +dx*(0:nbar-3))/90*( 14/3*(2+cosh(dx)) ...
                    + .5*(cosh(1.5*dx) + 9*cosh(1.25*dx) +23*cosh(.5*dx))...
                    +  1/6*(cosh(7/4*dx) + 121*cosh(.75*dx) +235*cosh(.25*dx)));

    Cons = 32*a^4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------     QUADRATIC  ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif order ==2  
    grand = @(w)rnCHF(w).*(sin(w/(2*a))./w).^3./(26*cos(w/a)+cos(2*w/a)+33);
    beta = real(fft([1/(960*a^3) exp(-1i*xmin*omega).*feval(grand,omega) ]));
    
    G = zeros(1,nbar+1);    
    G(nbar +1)  = W*(1/48 - exp(dx)*(exp(-11/8*dx)/720 + exp(-1.25*dx)/480 + exp(-9/8*dx)/80 + 7/1440*exp(-dx)));
    G(nbar)     = W*(.5 -.1*(7/24 + exp(-1.25*dx)/9 + exp(-dx)/6 + exp(-.75*dx) + 7/12*exp(-.5*dx) + 13/12*exp(-3/8*dx) + 11/24*exp(-dx/4) + 47/36*exp(-dx/8)));
    G(nbar-1)   = W*(47/48 - exp(-dx)*.1*(1 + exp(-1.25*dx)/9 + exp(-dx)/6 + exp(-.75*dx) + 7/9*exp(-.5*dx)...
                    + 44/9*cosh(dx/4) +7/12*exp(.5*dx) + 49/72*exp(5/8*dx) +3/16*exp(.75*dx) + 25/72*exp(7/8*dx) + 7/144*exp(dx)));
    G(1:nbar-2) = W - exp(xmin + dx*(0:nbar -3))*S_0*.1*(1 + 2/9*cosh(1.25*dx) +cosh(dx)/3 +2*cosh(.75*dx) + 14/9*cosh(.5*dx) +44/9*cosh(.25*dx));

    Cons = 960*a^(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------     LINEAR  ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif order == 1 
    grand = @(w)rnCHF(w).*(sin(w/(2*a))./w).^2./(2+cos(w/a));
    beta = real(fft([1/(24*a^2) exp(-1i*xmin*omega).*feval(grand,omega)]));
    
    G = zeros(1,nbar); 
    G(nbar)     = W*(.5 -(7/6 + 4/3*exp(-.75*dx) +exp(-.5*dx) + 4*exp(-.25*dx))/15);
    G(1:nbar-1) = W - exp(xmin + dx*(0:nbar-2))*S_0/15*(7/3 + 8/3*cosh(.75*dx) + 2*cosh(.5*dx)+8*cosh(.25*dx));
    
    Cons = 24*a^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------     HAAR  ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif order ==0 
    grand = @(w)rnCHF(w).*(sin(w/(2*a))./w);
    beta = real(fft([1/(4*a) exp(-1i*xmin*omega).*feval(grand,omega)]));

    G = zeros(1,nbar); 
    G(nbar)     = W*(.5 - a*(1-exp(-.5*dx)));
    G(1:nbar-1) = W- exp(xmin + dx*(0:nbar -2))*S_0*2*a*sinh(dx/2);
 
    Cons = 4*a;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if call == 1  % Use put-call parity
    price = Cons*exp(-r*T)/N*G*(beta(1:length(G))') + S_0*exp(-q*T) - W*exp(-r*T);
else
    price = Cons*exp(-r*T)/N*G*(beta(1:length(G))');
end

price = max(price, 0);  % Protect against deep out of money case


end

