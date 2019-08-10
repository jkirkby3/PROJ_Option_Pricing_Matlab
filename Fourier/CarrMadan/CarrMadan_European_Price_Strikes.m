function [ price ] = CarrMadan_European_Price_Strikes(S_0, W, rnCHF, N, T, r, q, call, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European Options using Method of Carr-Madan (1999)
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
% call  = 1 for call (else put)
% rnCHF = risk netural characteristic function (function handle with single argument)
%
% ----------------------
% Numerical Params 
% ----------------------
% N     = number of FFT grid points (power of 2, e.g. 2^12)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disc = exp(-r*T);

logK = log(W); % log of strikes vec
logS = log(S_0);
 
if nargin < 9
    alpha = .75;  % contour shift param (see Lee for recommendations)
end
eta = .05;

lam = 2*pi/(N*eta);  %spacing for log strike
b = N*lam/2; 

uv = 1:N;
ku = -b + lam*(uv - 1);
vj = (uv-1)*eta;
rnCHF = @(z) rnCHF(z).*exp(1i*z*logS);  % Convert to chf of log(S), rather than log return
Psij = feval(rnCHF, vj-(alpha+1)*1i)./(alpha^2 + alpha - vj.^2+1i*(2*alpha+1).*vj);

temp = (disc*eta/3)*exp(1i*vj*b).*Psij;
temp = temp.*(3+(-1).^uv - ((uv-1)==0));

Cku = real(exp(-alpha*ku).*fft(temp)/pi);

istrike = floor((logK + b)/lam +1);
xp = [ku(istrike) ku(istrike +1)];
yp = [Cku(istrike) Cku(istrike +1)];
price = real(interp1(xp, yp, logK));

if call ~=1  % If puts, price using Put-Call parity
    price = price - (S_0*exp(-q*T) - W*disc);
end

end

