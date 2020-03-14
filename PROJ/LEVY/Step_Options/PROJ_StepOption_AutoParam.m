function price = PROJ_StepOption_AutoParam(N,stepRho, call,down, S_0,W,H,M,r,q,rnCHF,T,L1,c2,c4, alphMult,TOLProb,TOLMean,rnCHFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for STEP-style barrier options and FADER options using PROJ method 
%     Step Payoff: exp(-stepRho*R) * (S_T - W)^+  for a call, where R is the proportion of time spent in knock-out region
%     Fader (Fade-in) Payoff: (1 - R) * (S_T - W)^+ for a call, where R is the proportion of time spent in knock-out region
%           (Fade-out) Can be priced by parity (Fade-in + Fade-out = Vanilla), so Price(Fade-out) = Price(Vanilla) - Price(Fade-in)
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% For algo/contract details and see section 4.4 of "Robust Barrier Option Pricing by Frame Projection under
%   Exponential Levy Dynamics", App. Math. Finance, 2017
% Contract is based on the paper "Step Options", (Linetsky, V.), Math. Finance 1999.
%
% ----------------------
% Contract/Model Params 
% ----------------------
% stepRho: if >= 0, then Step-Option: "softener" h(R) = exp(-stepRho * R)
%          if = -1, then Fader-Option: "softener" h(R) = 1 - R 
% S_0 = initial stock price (e.g. 100)
% W   = strike  (e.g. 100)
% r   = interest rate (e.g. 0.05)
% q   = dividend yield
% T   = time remaining until maturity (in years, e.g. T=1)
% M   = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
% call = 1 for call (else put)
% down = 1 for down and out (otherwise it's up and out)
% H    = barrier
% rnCHF, rnCHF_T = risk netural characteristic function (function handle with single argument),
%   at time steps dt=1/M and T, respectively
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% This version uses automated parameter selection, starting from the initial guess of N basis elements
% L1      = gridwidth param, e.g. L1 = 10 (grid width is iteratively increased if this is found to be insufficient)
% N       = initial choice of number of grid/basis points (power of 2, e.g. 2^12), auto increased if found to be insufficient
%           
% TOLProb = probability estimate accuracy, e.g. 5e-08;
% TOLMean = mean estimate accuracy, e.g. 1e-05;
% alphMult = used to increase grid width during parameter selection, e.g. alphMult = 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gamm = M+1;  % max number of time points (so max number that could possibly be spent in knockout region)
tauM = 1 / (M+1);

if stepRho >= 0 % Step Option  % NOTE: stepRho == 0 corresponds to vanilla option
    stepSoftener = @(u) exp(-stepRho*tauM*u);
    
elseif stepRho == -1  % Fader Option  (Fade-in)
    stepSoftener = @(u) 1 - u / (M + 1) ;
    
elseif stepRho == -2   % Ordinary Barrier option (no excursion forgiveness)
    stepSoftener = @(u) 1*(u==0);
end

gamm0 = 1; %HARDCODED: this param would allow us to specify an inital consumed budget

dt   = T/M;
nrdt = -r*dt;
h = log(H/S_0);
lws   = log(W/S_0);

%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;

alph = max(2*alphMult*abs(h), L1*sqrt(abs(c2*T) + sqrt(abs(c4*T))));


%=======================================
%Step 1: Satisfy Probability Tolerance
%=======================================
Nmax = max(2^8, 2*N);  %in case we already set N to be very large, we still need to enter both the loops

ErrProb = 10;
numTimesInLoop = 0;
N = N/2; alph = alph/alphMult;  %b/c we immediately double

while abs(ErrProb) > TOLProb && N < Nmax/2  %NOTE: strictly less here b/c we have to enter into the mean error loop at least once

    alph = alphMult*alph;
    N = 2*N;
    numTimesInLoop = numTimesInLoop + 1;
    
    fprintf('[-alpha,alpha] = [%.4f, %.4f] \n', -alph,alph)
    fprintf('\n')

    dx = 2*alph/(N-1); a = 1/dx;
    dw = 2*pi*a/N;
    xmin = -alph/2; xmax = xmin +(N/2 - 1)*dx;
    gam1 = (xmax - xmin)/2;  gam2 = (xmax + xmin)/2;
    grand  = dw*(1:(N-1));
    
    Prob = sum(exp(-1i*gam2*grand).*sin(gam1*grand)./grand.*rnCHFT(grand));
    Prob = Prob + sum(exp(1i*gam2*grand).*sin(gam1*grand)./grand.*rnCHFT(-grand)); %Note: sin(-x)/(-x) = sin(x)/x
    Prob = dw/pi*(gam1 + Prob);
    ErrProb = (1 - Prob);
    
    fprintf('NL = %.0f, ProbError: %.3e\n',numTimesInLoop,ErrProb)

end

%=======================================
%Step 2: Satisfy Mean Tolerance
%=======================================

numTimesInLoop = 0;
ErrMean = 10;
while abs(ErrMean) > TOLMean && N <= Nmax
    if numTimesInLoop>0  
        N = 2*N; %Only double N 
    end
    numTimesInLoop = numTimesInLoop + 1;
    
    dx = 2*alph/(N-1);

    xmin = -alph/2;    %ie log(S_0/S_0) - alph  ... note: gridwidth for log(S_T/S_0) is alph, so divide that by 2
    n_h = floor((h-xmin)/dx +1); 
    xmin = h - (n_h -1)*dx;    %realign so that h is on the grid (this is important for the case where h = 0)


    if h~= 0 %Realign so that h and 0 are both members of grid (if possible)
        nnot =  floor(1-xmin/dx);
        if abs(h) > dx  %so that n_h ~= nnot
            dx = (h - 0)/(n_h - nnot);
            xmin = dx*(1-nnot);  %hence nnot should remain on the grid
            %n_h = floor((h-xmin)/dx +1);  %NOT Numerically Stable
            n_h = floor(nnot + h/dx);  %Numerically Stable
        end
    else 
        nnot = n_h;  
    end

    a    = 1/dx;
    a2   = a^2;
    zmin = (1 - N/2)*dx;  %Kbar corresponds to zero
    
    Cons2 = 24*a2*exp(nrdt)/N;
    dw    = 2*pi*a/N;
    grand = dw*(1:N-1);
    grand = exp(-1i*zmin*grand).*rnCHF(grand).*(sin(grand/(2*a))./grand).^2./(2+cos(grand/a));
    beta  = Cons2*real(fft([1/(24*a2) grand]));   %%%%  NOTE: all toep matrices incorporate exp(-r*dt)


    varthet_star = (exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4) + exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4))/18;

    ErrMean = exp(r*dt)*varthet_star*(beta*exp(zmin + dx*(0:N-1))') - exp((r-q)*dt);    %note: multiplying by exp(r*dt) to cancel the exp(-r*dt) in Cons2
    ErrMean = ErrMean*M*S_0;
    fprintf('NLMean = %.0f, ErrMean = %.3e\n',numTimesInLoop, ErrMean)
end
fprintf('Final log2(N): %.0f \n',log2(N))

interp_Atend = 0;
if 0 < abs(h) && abs(h)<dx
    interp_Atend = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   DETERMINE COMMON Params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K     = N/2;
nbar  = floor(a*(lws - xmin)+1);
rho   = lws - (xmin+(nbar - 1)*dx);
zeta  = a*rho;

toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];   toepM = fft(toepM);

%%%% PAYOFF CONSTANTS-----------------------------------
varthet_01 = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
varthet_m10 = exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4)/18;
varthet_star = varthet_01 + varthet_m10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if down == 1 && call ~= 1  %DOP
    %l = log(H/S_0);
    n_l = n_h;

    zeta_plus = zeta*q_plus; zeta_minus = zeta*q_minus;
    rho_plus = rho*q_plus; rho_minus = rho*q_minus;

    ed1 = exp(rho_minus); ed2 = exp(rho/2); ed3 = exp(rho_plus);

    dbar_1 = zeta^2/2;
    dbar_0 = zeta - dbar_1;         %  dbar_1 = zeta + .5*((zeta - 1)^2 - 1);
    d_0    = zeta*(5*( (1-zeta_minus)*ed1 + (1-zeta_plus)*ed3 ) + 4*(2-zeta)*ed2)/18;
    d_1    = zeta*( 5*(zeta_minus*ed1 + zeta_plus*ed3) + 4*zeta*ed2 )/18;            

    Thet = zeros(K,1);
    Thet(1:nbar-1) =  W - exp(xmin +dx*(0:nbar-2))*S_0*varthet_star;
    Thet(nbar)     =  W*(.5 + dbar_0 - exp(-rho)*(varthet_m10 + d_0));
    Thet(nbar + 1) =  W*(dbar_1 - exp(- rho)*d_1);
  
    %%%%%%% Initialize Val
    Val = zeros(K,Gamm+1);  % j = 1 corresponds to gamm = 0

    ThetTemp = zeros(K,1);
    for j=1:Gamm
       gamm = j-1;
       ThetTemp(1:n_l-1) = stepSoftener(gamm+1)*Thet(1:n_l-1);
       
       ThetTemp(n_l) = .5*(stepSoftener(gamm+1) + stepSoftener(gamm))*Thet(n_l);
       
       ThetTemp(n_l+1:K) = stepSoftener(gamm)*Thet(n_l+1:K);
       
       p   = ifft(toepM.*fft([ThetTemp;zeros(K,1)]));
       Val(:,j) = p(1:K); 
    end

    %  NOTE: this introduces a kink and we use interpolation in next step
    %  ... this can be improved 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% CUMULATIVE PARISIAN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %NOTE: DIFFERS FROM parisian in the first step ... M-1
    for m=M-2:-1:0

        for j=1:m+1 %NOTE: dont need to compute all Gamma of them
            Thet(1)      = (13*Val(1,j+1)+15*Val(2,j+1)-5*Val(3,j+1)+Val(4,j+1))/48;
            Thet(2:n_l-1) = (Val(1:n_l-2,j+1)+10*Val(2:n_l-1,j+1)+Val(3:n_l,j+1))/12;

            Thet(n_l)     = (13*Val(n_l,j+1)+15*Val(n_l-1,j+1)-5*Val(n_l-2,j+1)+Val(n_l-3,j+1))/48 ...
                          + (13*Val(n_l,j)+15*Val(n_l+1,j)-5*Val(n_l+2,j)+Val(n_l+3,j))/48;                          

            Thet(n_l+1:K-1) = (Val(n_l:K-2,j)+10*Val(n_l+1:K-1,j)+Val(n_l+2:K,j))/12;
            Thet(K)        = (13*Val(K,j)+15*Val(K-1,j)-5*Val(K-2,j)+Val(K-3,j))/48;

            p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
            Val(:,j)  = p(1:K);
        end

        %%% Now to Gamm+1
        j = Gamm+1;
        Thet(1:n_l-1)   = 0;
        Thet(n_l)       = (13*Val(n_l,j)+15*Val(n_l+1,j)-5*Val(n_l+2,j)+Val(n_l+3,j))/48;
        Thet(n_l+1:K-1) = (Val(n_l:K-2,j)+10*Val(n_l+1:K-1,j)+Val(n_l+2:K,j))/12;

        p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
        Val(:,j)  = p(1:K);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif down ~=1 && call == 1  %UOC
    %u    = log(H/S_0);
    n_u = n_h;

    sigma = 1 - zeta; sigma_plus = (q_plus-.5)*sigma; sigma_minus = (q_minus-.5)*sigma;

    es1 = exp(dx*sigma_plus); es2 = exp(dx*sigma_minus);

    dbar_0 = .5 + zeta*(.5*zeta-1);
    dbar_1 = sigma*(1 - .5*sigma);

    d_0 = exp((rho+dx)*.5)*sigma^2/18*(5*((1-q_minus)*es2 +(1-q_plus)*es1) + 4);
    d_1 = exp((rho-dx)*.5)*sigma/18*(5*( (.5*(zeta+1) +sigma_minus)*es2 + (.5*(zeta+1) +sigma_plus)*es1 ) + 4*(zeta+1) );
           
    Thet = zeros(K,1);
    Thet(nbar)     =  W*(exp(-rho)*d_0 - dbar_0);
    Thet(nbar + 1) =  W*(exp(dx-rho)*(varthet_01 +d_1) -(.5 + dbar_1) );
    Thet(nbar + 2:K)  = exp(xmin +dx*(nbar+1:K-1))*S_0*varthet_star - W;

    Val = zeros(K,Gamm+1);
    
    ThetTemp = zeros(K,1);
    for j=1:Gamm
       gamm = j-1;
       ThetTemp(1:n_u-1) = stepSoftener(gamm)*Thet(1:n_u-1);
       
       ThetTemp(n_u) = .5*(stepSoftener(gamm) + stepSoftener(gamm+1))*Thet(n_u);
       
       ThetTemp(n_u+1:K) = stepSoftener(gamm+1)*Thet(n_u+1:K);
       
       p   = ifft(toepM.*fft([ThetTemp;zeros(K,1)]));
       Val(:,j) = p(1:K); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% CUMULATIVE PARISIAN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for m=M-2:-1:0
        for j=1:Gamm
            Thet(1)      = (13*Val(1,j)+15*Val(2,j)-5*Val(3,j)+Val(4,j))/48;
            Thet(2:n_u-1) = (Val(1:n_u-2,j)+10*Val(2:n_u-1,j)+Val(3:n_u,j))/12;

            Thet(n_u)     = (13*Val(n_u,j)+15*Val(n_u-1,j)-5*Val(n_u-2,j)+Val(n_u-3,j))/48 ...
                          + (13*Val(n_u,j+1)+15*Val(n_u+1,j+1)-5*Val(n_u+2,j+1)+Val(n_u+3,j+1))/48;

            Thet(n_u+1:K-1) = (Val(n_u:K-2,j+1)+10*Val(n_u+1:K-1,j+1)+Val(n_u+2:K,j+1))/12;
            Thet(K)        = (13*Val(K,j+1)+15*Val(K-1,j+1)-5*Val(K-2,j+1)+Val(K-3,j+1))/48;

            p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
            Val(:,j)  = p(1:K);
        end

        %%% Now to Gamm+1
        j = Gamm+1;
        Thet(1)       = (13*Val(1,j)+15*Val(2,j)-5*Val(3,j)+Val(4,j))/48;
        Thet(2:n_u-1) = (Val(1:n_u-2,j)+10*Val(2:n_u-1,j)+Val(3:n_u,j))/12;
        Thet(n_u)     = (13*Val(n_u,j)+15*Val(n_u-1,j)-5*Val(n_u-2,j)+Val(n_u-3,j))/48;
        Thet(n_u+1:K)  =0;

        p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
        Val(:,j)  = p(1:K);
    end
end


if interp_Atend ~= 1
    price = Val(nnot,gamm0);

else  %%% INTERPOLATION
    dd = 0 - (xmin+ (nnot -1)*dx); 
    price = Val(nnot,gamm0) + (Val(nnot+1,gamm0) - Val(nnot,gamm0))*dd/dx;
end


end

