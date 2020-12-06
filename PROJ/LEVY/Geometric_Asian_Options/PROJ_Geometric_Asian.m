function price = PROJ_Geometric_Asian(N, alph, S_0, M, W, call, T, r, q, rnSYMB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Geometric Asian Options using PROJ method
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% Reference:   (1) An Efficient Transform Method For Asian Option Pricing, SIAM J. Financial Math., 2016
%              (2) Efficient Option Pricing by Frame Duality with the Fast Fourier Transform. 
%                  SIAM J. Financial Math (2015), Kirkby, J.L
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
% rnSYMB =  risk neutral symbol of log return over time step dt = 1/M (function handle with single argument)
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% alph  = grid with is 2*alph
% N     = number of grid/basis points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Psi = @(u) rnSYMB(u) - (r-q-rnSYMB(-1i))*u;  % Un-risk neutralize it

dx = 2*alph/(N-1); a = 1/dx;
dt    = T/M;

xmin  = log(W)-dx;
if abs(xmin) > 1.01*alph
    alph = 1.01*abs(xmin);
    price = PROJ_Geometric_Asian(N, alph, S_0, M, W, call, T, r, q, rnSYMB);
    return;
end

dw    = 2*pi*a/N;
omega = (dw: dw: (N-1)*dw);  %We calcuate coefficient of w=0 explicitly
gran  = zeros(1,N-1);

for j=1:(N-1)
    gran(j) = sum(Psi(omega(j)*(1-(1:M)/(M+1))));
end

a2 = a^2; a3 = a*a2;

gg   = @(w)(sin(w/(2*a))./w).^3./(26*cos(w/a)+cos(2*w/a)+33);
gran = exp(1i*(log(S_0) +.5*T*(r-q-Psi(-1i)))*omega + dt*gran).*gg(omega);
beta = real(fft([1/(960*a3) exp(-1i*xmin*omega).*gran ]));

x1   = xmin + dx;   
ex1  = exp(x1); 
G    = zeros(1,N/2);
G(1) = ex1*(a3*exp(.5*dx)-a*(1/8+a/2+a2)) - W/48;
G(2) = ex1*(a3*(2+exp(1.5*dx)-3*exp(.5*dx))-.75*a) - .5*W;
G(3) = ex1*(.5*a2 -a/8 +a3*(exp(2.5*dx)+3*(exp(.5*dx)-exp(1.5*dx))-1)) - 47/48*W;
G(4:N/2) = exp(x1+dx*(2:N/2-2))*a3*(2*sinh(1.5*dx)-6*sinh(.5*dx))-W;


price = 960*a^(3)*exp(-r*T)/N*G*beta(1:N/2)';
price = max(0, price);
if ~call==1
    error("Sorry, havent yet added the put option... just use PCP");
end


end

