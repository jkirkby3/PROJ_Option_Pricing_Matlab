function price = Bermudan_PUT_func_alpha(M, S_0, W, r, T, rnCHF, N, alph)
%UNTITLED8 Summary of this function goes here
%   Note: q enters through rnCHF, but is not otherwise needed
% NOTE: uses no augmentation

dt   = T/M;

K = N/2;

dx = 2*alph/(N-1); a = 1/dx;

ThetM = zeros(K,1); 
Cons3 = 1/48;
Cons4 = 1/12;   

lws = log(W/S_0); % Center value grid on [lws - 2^(Pbar-1),lws + 2^(Pbar-1)]
% Will perturb so that log(S_0/S_0)=0 is a member
% note: xnbar = lws


%%% TRY NEW GRID APPROACH 12/29/2015
nnot = K/2;
dxtil = 1/a; %for now... change to dxtil = 2*alpha/(N-1)
nbar = floor(lws*a +K/2);
if abs(lws)<dxtil
    dx = dxtil;  
elseif lws<0
    dx = lws/(1+nbar - K/2);
    nbar = nbar+1;
elseif lws>0
    dx = lws/(nbar-K/2);

end
a = 1/dx;
xmin = (1-K/2)*dx;


a2    = a^2;  
Cons2 = 24*a2*exp(-r*dt)/N;
zmin  = (1 - K)*dx;  %Kbar corresponds to zero

dw    = 2*pi*a/N;
grand = (dw: dw: (N-1)*dw);
grand = exp(-1i*zmin*grand).*rnCHF(grand).*(sin(grand/(2*a))./grand).^2./(2+cos(grand/a));
beta  = Cons2*real(fft([1/(24*a2) grand]));   %%%%  NOTE: all toep matrices incorporate exp(-r*dt)

toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];
toepM = fft(toepM);

%%%%%%%
Gs  = zeros(K,1);
Gs(1:nbar) = exp(xmin + dx*(0:nbar-1))*S_0;


%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;

%%%% PAYOFF CONSTANTS-----------------------------------
varthet_01 = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
varthet_m10 = exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4)/18;
varthet_star = varthet_01 + varthet_m10;

ThetM(nbar)     = W*(.5 - varthet_m10);
ThetM(1:nbar-1) = W - varthet_star*Gs(1:nbar-1);
Gs(1:nbar) = W - Gs(1:nbar);

% %%%% Compute Augmentation
% toepL = [zeros(K,1); 0 ; beta(K-1:-1:1)'];
% toepL = fft(toepL);
% Thetbar2 = W - S_0*varthet_star*exp(xmin - dx*(K:-1:1))';
% p = ifft(toepL.*fft([Thetbar2; zeros(K,1)]));
% theta_aug  = p(1:K);  %already includes the exp(-r*dt) through Beta

%%%%%%%
p = ifft(toepM.*fft([ThetM(1:K);zeros(K,1)]));
Cont = p(1:K); %+ theta_aug;
%%%%%%%

Thet = zeros(K,1);

kstr = nbar +1;
    
for m=M-2:-1:0
    
    while kstr>1 && Cont(kstr)>Gs(kstr)
        kstr = kstr -1;
    end
    if kstr>=2
        xkstr = xmin +(kstr -1)*dx;

        Ck1 = Cont(kstr-1); Ck2 = Cont(kstr); Ck3 = Cont(kstr+1);
        %%% Linear interp of payoff
        Gk2 = Gs(kstr); Gk3 = Gs(kstr+1);
        %%% xstr
        tmp1 = Ck2-Gk2; tmp2 = Ck3 - Gk3;
        xstrs  = ((xkstr+dx)*tmp1 - xkstr*tmp2)/(tmp1-tmp2);
    else
        xkstr = xmin;
        kstr = 1; xstrs = xmin;
        Ck2 = Cont(kstr); Ck1 = Ck2;  Ck3 = Cont(kstr+1);
    end
        
    rho = xstrs-xkstr;
    zeta = a*rho;
    
    zeta2 = zeta^2; zeta3 = zeta*zeta2; zeta4 = zeta*zeta3;
    
    zeta_plus = zeta*q_plus; zeta_minus = zeta*q_minus;
    rho_plus = rho*q_plus; rho_minus = rho*q_minus;

    ed1 = exp(rho_minus); ed2 = exp(rho/2); ed3 = exp(rho_plus);

    dbar_1 = zeta2/2;
    dbar_0 = zeta - dbar_1;         %  dbar_1 = zeta + .5*((zeta - 1)^2 - 1);
    d_0    = zeta*(5*( (1-zeta_minus)*ed1 + (1-zeta_plus)*ed3 ) + 4*(2-zeta)*ed2)/18;
    d_1    =  exp(-dx)*zeta*( 5*(zeta_minus*ed1 + zeta_plus*ed3) + 4*zeta*ed2 )/18;    
    
    
    Thet(1:kstr -1) = ThetM(1:kstr-1);
    
    Ck4 = Cont(kstr +2);
    Thet(kstr) = W*(.5 + dbar_0) -S_0*exp(xkstr)*(varthet_m10 +d_0)...
                 + zeta4/8*(Ck1 -2*Ck2 + Ck3) +zeta3/3*(Ck2-Ck1)...
                 + zeta2/4*(Ck1 +2*Ck2 - Ck3) - zeta*Ck2...
                 -Ck1/24 +5/12*Ck2 +Ck3/8;
     
    Thet(kstr+1) = W*dbar_1 - S_0*exp(xkstr+dx)*d_1    + zeta4/8*(-Ck2 +2*Ck3 - Ck4)...
                  + zeta3/6*(3*Ck2 - 4*Ck3 + Ck4) - .5*zeta2*Ck2...
                  + Cons4*(Ck2 +10*Ck3 + Ck4);
    Thet(kstr+2:K-1) = Cons4*(Cont(kstr+1:K-2)+10*Cont(kstr+2:K-1)+Cont(kstr+3:K));
    Thet(K)         = Cons3*(13*Cont(K)+15*Cont(K-1)-5*Cont(K-2)+Cont(K-3));
    
    p         = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));

    Cont(1:K) = p(1:K); %+ theta_aug;
end

price = Cont(nnot);

end

