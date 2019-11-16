function Value = PROJ_Lookback( N, alpha, S_0, W, call, r, q, M, T, rnSYMB, floating_strike)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Discrete Lookback/Hindsight Options using PROJ method
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
% floating_strike = 1 for floating strike contract (else fixed strike)
% rnSYMB = risk netural symbol function (function handle with single argument)
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% alpha  = grid with is 2*alpha
% N     = number of grid/basis points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = 2*alpha/(N-1); a = 1/dx;
K    = N/2;  
dt   = T/M;
zmin  = (1 - K)*dx; 

if (floating_strike == 1 && call ~= 1) || (floating_strike ~= 1 && call == 1)
    % floating strike put or fixed strike call (ie we find max of process)
    rnCHFstar = @(u)exp(dt*(rnSYMB(u-1i) - (r-q))); 
else
    error("This contract type not currently finished. Implementaion in progress");
    %rnCHFstar = @(u)exp(dt*(rnSYMB(-u+1i) + (r-q)));   % min(X) = -max(-X)
    %rnCHFstar = @(u)exp(dt*(rnSYMB(u-1i) + (r-q)));   % min(X) = -max(-X)
    
    rnCHFstar = @(u)exp(dt*(rnSYMB(-u-1i) - (r-q))); 
end

Cons   = 32*a^4/N;   
b0     = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
grandG = @(w)rnCHFstar(w).*(sin(w/(2*a))./w).^4./(b0 + b1*cos(w/a) +b2*cos(2*w/a) +b3*cos(3*w/a));
dw    = 2*pi*a/N; grand = dw*(1:N-1);   %dw     = 2*pi*a/N; grand = dw*(1:N-1);
beta   = Cons*real(fft([1/(32*a^4) exp(-1i*zmin*grand).*feval(grandG,grand)]));

Lambda = a*(2/3*beta(K:-1:1)' + 1/6*(beta(K+1:-1:2)' + [beta(K-1:-1:1) 0]'));  
toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];  toepM = fft(toepM);

%%%% Simpson rule
weights = ones(1,K); weights(1) = .5; weights(K)=.5;

%%%% m = 1
fvec = Lambda;  %initialize fvec
Thet = zeros(K,1);

%%%% Extra
Evec1 = [beta(K-1:-1:1) 0 ]';

d11 = 1757/17280;  d22 = 21757/30240; d33 = 2371/20160;  d44 = 149/6048;  d55 = -613/120960;
h2 = -1/720; h1 = 31/180; h0 = 79/120;

%Cubic interpolation for the extra vec
v1 = 73/2520; v2 = 47/2520; v3 = -19/2520; v4 = 1/630; 

for m=2:1:M
    zm = 1 - dx*weights*fvec;   %%%This is really z_{m-1}
    
    %%%%%%%%%%%%
    Thetastar1  = v1*fvec(1) + v2*fvec(2) + v3*fvec(3) + v4*fvec(4);
    Thet(1)     = (13*fvec(1)+15*fvec(2)-5*fvec(3)+fvec(4))/48;
    Thet(2)     = d11*fvec(1) + d22*fvec(2) + d33*fvec(3) + d44*fvec(4) + d55*fvec(5);
    Thet(3:K-2) = h0*fvec(3:K-2) + h1*(fvec(2:K-3)+fvec(4:K-1)) +h2*(fvec(1:K-4)+fvec(5:K));
    Thet(K-1)   = d11*fvec(K) + d22*fvec(K-1) + d33*fvec(K-2) + d44*fvec(K-3) + d55*fvec(K-4);
    Thet(K)     = (13*fvec(K)+15*fvec(K-1)-5*fvec(K-2)+fvec(K-3))/48;
    %%%%%%%%%%%%
    
    p = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
    fvec(1:K) = p(1:K)+ zm*Lambda + Thetastar1*Evec1;
end


% All contracts are priced w.r.t floating strike looback put using parity

if floating_strike == 1
    if call ~= 1  % Floating strike Lookback put:  max{S_m: 0<=m<=M} - S_T
        %===========================
        E = (dx*(weights*( (exp(dx*(0:K-1)')-ones(K,1)).*fvec)));   % this is: E[exp(Y_M)] - 1
        Value = S_0*exp(-q*T)*(E);   % V = S_0*exp(-r*T) E[exp(M)] - exp(-q*T)*S_0 = S_0*exp(-q*T)*(E[exp(Y_M)] - 1)
        %===========================
        
    else  % Floating strike Lookback call: S_T - min{S_m: 0<=m<=M}
        E = (dx*(weights*( (ones(K,1) - exp(-dx*(0:K-1)')).*fvec))); 
        Value = min(W,S_0) *exp(-q*T)*(E);
        
    end
else % fixed strike (hindsight)
    
    if call ~=1 % Fixed Strike (Hindsight) Put: (W - min{S_m: 0<=m<=M})^+
        % W*exp(-r*T) - S_0*exp(-q*T) + 
        E = (dx*(weights*( (ones(K,1) - exp(-dx*(0:K-1)')).*fvec))); 
        V_fsc = min(W,S_0) *exp(-q*T)*(E);
        Value = W*exp(-r*T) - S_0*exp(-q*T) + V_fsc;
        
    else  % Fixed Strike (Hindsight) Call: (max{S_m: 0<=m<=M} - W)^+
        %===========================
        E = (dx*(weights*( (exp(dx*(0:K-1)')-ones(K,1)).*fvec)));   % this is: E[exp(Y_M)] - 1
        V_fsp = max(S_0, W)*exp(-q*T)*(E);   % floating strike put with different starting value
        Value = V_fsp + S_0*exp(-q*T) - exp(-r*T)*W;
        %===========================
    end
end


end