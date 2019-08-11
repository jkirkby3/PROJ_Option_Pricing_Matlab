function [ price ] = CONV_European_Price(S_0, W, rnCHF, T, r, call, N, alph, damp_alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for European Options using the CONV method of Lord, Fang, Bervoets, and Oosterlee (2008)
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
% N  = number of FFT grid points (power of 2, e.g. 2^12)
% alpha = gridwidth param, density centered on [-alph, alph], determined e.g. using cumulants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if call ~= 1
    call = -1;
end

if nargin < 9
    damp_alpha = -call*0.5;   % damping factor
end
                                    
dx = (2*alph)/N;             %  density centered on [-alph, alph]                                                
du = 2 * pi / (2*alph);                           

grid1 = (0:N-1)';                          
grid2 = (-1).^grid1;   

x = (grid1 .* dx) - (N/2 * dx); 
y = log(W / S_0) + (grid1 .* dx) - (N/2 * dx); 

u = (grid1 .* du) - (N/2 * du); 

payoff = max(call*(S_0*exp(y) - W), 0);    
damped = payoff .* exp(damp_alpha .* y);   
w = ones(N,1); w(1) = 0.5; w(N) = 0.5;  % trapezoidal quadrature
f = ifft( (grid2) .* w .* damped );    

chfval = rnCHF(-(u - (1i*damp_alpha)));

f = exp( 1i .* grid1 .* (y(1) - x(1)) .* du ).* chfval.* f;    

% Final Valuation Step
C = abs(exp( -(r * T) - (damp_alpha .* x) + (1i .* u .* (y(1) - x(1))) ).* grid2 .* fft(f));        
price = double(C(N/2 + 1, 1));       


end

