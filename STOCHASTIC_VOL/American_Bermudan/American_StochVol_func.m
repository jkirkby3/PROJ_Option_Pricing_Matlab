function price = American_StochVol_func( N,alph,M,r,T,S_0,W,m_0,psi_J,model, modparam, gridMethod, gamma, gridMultParam)
%  American PUT Option
%-----------------------------
% N  : size of density grid (value grid is K:=N/2)
% alph: density gridwith param, density on [-\alph,\alph]
% r  : interest rate (note, can also add dividend yield .. needs to be part of drift below, ie r-q)
% T  : number of years (T = 2 is two years, T = .5 is half a year)
% S_0: initial Underlying
% M  : number of subintervals of [0,T] (total of M+1 points in time grid)
% W  : strike  (used instead of K)
% m_0: number of states to approximate the Heston model with
% psi_J: characteristic exponenent of jump part...
%        function handdle: psi_J(xi) = lambda*(phi(xi) -1)
% model: 1 = Heston, 2 = SteinStein, 3 = 3/2 Model, 4 = 4/2 Model, 
%        5 = HullWhite, 6 = Scott, 7 = alpha-Hypergeometric
% gamma: var grid width parameter, grid is +/- gamma*stddev(variance process)
% gridMethod: 1 for uniform, 4 for nonuniform (the sinh approach.. assumes v0 as center)
% modparam: contains all necessary params for the specific model (see below during assingment which ones are needed)
%-------------------------------
%%% Note: be careful about the parameter Rho (vs rho used in algorithm)

K    = N/2;
dx   = 2*alph/(N-1);
lws  = log(W/S_0);
dt   = T/M;

%%% GRID which aligns 0 as well as log(W/S_0)
nnot = K/2;
dxtil = dx; %for now... change to dxtil = 2*alpha/(N-1)
nbar = floor(lws/dx +K/2);
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


%%%%////////////////////////////////////////////////////////
%%%% Initialize THETA ... 
%%%%////////////////////////////////////////////////////////
THET  = zeros(K,m_0);  %Each column corresponds to a regime state

%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;


%%%% PAYOFF CONSTANTS-----------------------------------
varthet_01 = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
varthet_m10 = exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4)/18;
varthet_star = varthet_01 + varthet_m10;

%%%%%%%
Gs  = zeros(K,1);
Gs(1:nbar) = exp(xmin + dx*(0:nbar-1))*S_0;

%Define Terminal Theta Coeffs 
ThetM = zeros(K,1);
ThetM(1:nbar-1) = W - varthet_star*Gs(1:nbar-1);
ThetM(nbar)     = W*(.5 - varthet_m10);
Gs(1:nbar)      = W - Gs(1:nbar); %For American options (otherwise no need for Gs to be defined separately)


%%%%////////////////////////////////////////////////////////
%%%% Intialize Q matrix and variance set
%%%%////////////////////////////////////////////////////////
t = T/2;

if model == 1 %HESTON  (eta, theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; 
    mu_func  = @(v)eta*(theta-v);  %variance process drift
    sig_func = @(v)Sigmav*sqrt(v); %variance process vol
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmav^2/eta*v0*(exp(-eta*t)-exp(-2*eta*t)) +theta*Sigmav^2/(2*eta)*(1-exp(-eta*t)+exp(-2*eta*t));
    
elseif model == 2 || model == 6 %STEIN STEIN / SCOTT (eta, theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0;
    mu_func  = @(v)eta*(theta-v);
    sig_func = @(v)Sigmav*(v>-100); %note the (v>-100) just ensures that sig_func returns a vector
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmav^2/(2*eta)*(1 -exp(-2*eta*t));
    
elseif model == 3 % 3/2 MODEL
    %Transform to parameters that can be use in 4/2 model
    eta = modparam.eta*modparam.theta;
    theta = (modparam.eta + modparam.Sigmav^2)/eta;
    Sigmav = -modparam.Sigmav;
    v0 = 1/modparam.v0;
    Rho = modparam.rho;
    aa = 0; bb = 1; %Now use the 4/2 Model:
    
    mu_func  = @(v)eta*(theta-v);
    sig_func = @(v)Sigmav*sqrt(v);
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmav^2/eta*v0*(exp(-eta*t)-exp(-2*eta*t)) +theta*Sigmav^2/(2*eta)*(1-exp(-eta*t)+exp(-2*eta*t));
    
elseif model == 4 % 4/2 MODEL (aa,bb, eta,theta, Rho, Sigmav, v0)
    eta = modparam.eta; theta = modparam.theta; Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0;
    aa = modparam.aa; bb = modparam.bb;
    
    mu_func  = @(v)eta*(theta-v);
    sig_func = @(v)Sigmav*sqrt(v);
    
    mu_H = exp(-eta*t)*v0 + theta*(1-exp(-eta*t));
    sig2_H = Sigmav^2/eta*v0*(exp(-eta*t)-exp(-2*eta*t)) +theta*Sigmav^2/(2*eta)*(1-exp(-eta*t)+exp(-2*eta*t));
    
elseif model == 5 % Hull White Model  (av, Rho, Sigmav, v0)
    Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; av = modparam.av;
    
    mu_func  = @(v)av*v;
    sig_func = @(v)Sigmav*v;
    
    mu_H = v0*exp(av*t);
    sig2_H = v0^2*exp(2*av*t)*(exp(Sigmav^2*t)-1);
elseif model == 7 % alpha-Hypergeometric
%      %%% Commented Version is With the Second Formultion
%     Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; alphav = modparam.alphav; av = modparam.av; bv = modparam.bv;
%     
%     mu_func = @(v)2*(a+Sigmav^2)*v - 2*bv*v.^(1+alphav/2);
%     sig_func = @(v)2*Sigmav*v;
%     
%     Av     = 2*(av + Sigmav^2 - (v0)^(alphav/2)); 
%     mu_H   = v0*exp(Av*t);
%     sig2_H = v0^2*exp(2*av*t)*(exp((2*Sigmav)^2*t)-1);
    Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; eta = modparam.eta; av = modparam.av; theta = modparam.theta;
    mu_func = @(v)eta - theta*exp(av*v);
    sig_func = @(v)Sigmav*(v>-100);
    
    %%% Estimate Mean and Variance Using Stein Stein and First order approx
    EtaBar = theta*av*exp(av*v0); ThetaBar = (eta - theta*exp(av*v0)*(1-av*v0))/EtaBar;
    mu_H = exp(-EtaBar*t)*v0 + ThetaBar*(1-exp(-EtaBar*t));
    %mu_H = v0;
    sig2_H = Sigmav^2/(2*EtaBar)*(1 -exp(-2*EtaBar*t));
    
elseif model == 8 % JACOBI
    Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; eta = modparam.eta; theta = modparam.theta; vmin = modparam.vmin; vmax = modparam.vmax;
    %Qfunc = @(v) (v - vmin).*(vmax - v)/denomQ;
    %Qsqrt = @(v) sqrt((v - vmin).*(vmax - v)/denomQ);
    denomQ = (sqrt(vmax) - sqrt(vmin))^2;
    mu_func  = @(u) eta*(theta - u);
    sig_func = @(u) Sigmav/sqrt(denomQ)*sqrt((u - vmin).*(vmax - u));
end

%%% Determine Variance Grids
if model == 8 %JACOBI
    lx = vmin ; %For now (but maybe make it a little bigger than vmin, so vol doesnt hit zero)
    ux = vmax; %For now
else
    if model == 2 %For Stein-Stein, more sensitive
        lx = max(0.01,mu_H - gamma*sqrt(sig2_H));  %variance grid lower bound
    elseif model == 6 || model == 7   %For Scott/alpha-Hyper, allow lx negative
        lx = mu_H - gamma*sqrt(sig2_H);
    else
        lx = max(0.00001,mu_H - gamma*sqrt(sig2_H));  %variance grid lower bound
    end
    ux = mu_H + gamma*sqrt(sig2_H);  %variance grid upper bound
end

boundaryMethod = 1;
center = v0; %this is where grid clusters... we can experiment with other choices.. 
[Q,v]  = Q_Matrix_AllForms(m_0,mu_func,sig_func,lx,ux,gridMethod, gridMultParam, center, boundaryMethod);


% GridMultParam = 0.2;
% [Q,v]  = General_Q_Matrix_Newest(m_0,mu_func,sig_func,lx,ux,gridMethod,center, GridMultParam);

%%%%////////////////////////////////////////////////////////
%%%% Populate the Matrix Exponentials
%%%%////////////////////////////////////////////////////////

dxi    = 2*pi*a/N;
xi     = dxi*(0:N-1)';  

if model == 1 %HESTON
    c1 = (Rho*eta/Sigmav - .5);  c2 = (r - Rho*eta*theta/Sigmav);   c3 = .5*(1-Rho^2);
    v1 = dt*1i*(c1*v + c2 - psi_J(-1i));  %Note: we now have the compensated jump component
    v2 = dt*c3*v;
    fv = (1i*dxi*Rho/Sigmav)*v;
    
elseif model == 2 %STEIN STEIN
    c1 = Rho*eta/Sigmav - .5;  c2 = Rho*eta*theta/Sigmav;  c3 = r  - psi_J(-1i) - Rho*Sigmav/2; 
    vsq = v.^2;
    v1 = dt*1i*(c1*vsq - c2*v + c3);
    v2 = dt*.5*vsq*(1-Rho^2);
    fv = (1i*dxi*.5*Rho/Sigmav)*vsq;
    
elseif model == 3 || model  == 4 % 3/2 or 4/2 Model  (3/2 uses transformation then applies 4/2)
    c1 = aa*Rho*eta/Sigmav - .5*aa^2; c2 = .5*(Rho*bb*Sigmav-bb^2) - bb*Rho*eta*theta/Sigmav;
    c3 = Rho*eta/Sigmav*(bb-aa*theta) + r - psi_J(-1i) - aa*bb; 
    sqrtv = sqrt(v);
    v1 = dt*1i*(c1*v + c2./v + c3);
    v2 = dt*.5*(1-Rho^2)*(aa*sqrtv + bb./sqrtv).^2;
    fv = (1i*dxi*Rho/Sigmav)*(aa*v + bb*log(v));

elseif model == 5 % Hull White Model
    c1 = .25*Rho*Sigmav - av*Rho/Sigmav;  c2 = .5;  c3 = r - psi_J(-1i);
    sqrtv = sqrt(v);
    v1 = dt*1i*(c1*sqrtv - c2*v + c3);
    v2 = dt*.5*(1-Rho^2)*v;
    fv = 1i*dxi*2*Rho/Sigmav*sqrtv;
elseif model == 6 % Scott
    c1v = Rho*(eta/Sigmav*(v - theta)-Sigmav/2);  %Note: this depends on v
    c2 = .5; c3 = r - psi_J(-1i);
    expv = exp(v); expv2 = expv.^2;
    v1 = dt*1i*(c1v.*expv - c2*expv2 + c3);
    v2 = dt*.5*(1-Rho^2)*expv2;
    fv = 1i*dxi*Rho/Sigmav*expv;
elseif model == 7 % alpha-Hypergeometric
    %%% Commented Version is With the Second Formultion
%     c1 = (Rho*Sigmav*.5 - Rho*(av + Sigmav)/Sigmav);
%     c2 = .5; c3 = Rho*bv/Sigmav; c4 = r - psi_J(-1i);
%     sqrtv = sqrt(v);
%     v1 = dt*1i*(c1*sqrtv  - c2*v + c3*v.^(1+alphav)/2 + c4 );
%     v2 = dt*.5*(1-Rho^2)*v;
%     fv = 1i*dxi*Rho/Sigmav*sqrtv;

    c1 = Rho*theta/Sigmav;
    c2 = Rho*(eta/Sigmav + Sigmav/2); c3 = .5; c4 = r - psi_J(-1i);
    expv = exp(v); expv2 = expv.^2;
    v1 = dt*1i*(c1*exp((1+av)*v) - c2*expv - c3*expv2 + c4);
    v2 = dt*.5*(1-Rho^2)*expv2;
    fv = 1i*dxi*Rho/Sigmav*expv;
elseif model == 8 % JACOBI
%     c1 = (Rho*eta/Sigmav - .5);  c2 = (r - Rho*eta*theta/Sigmav);   c3 = .5*(1-Rho^2);
%     v1 = dt*1i*(c1*v + c2 - psi_J(-1i));  %Note: we now have the compensated jump component
%     v2 = dt*c3*v;
    c1 = r - Rho*eta*theta/Sigmav - psi_J(-1i);
    c2 = Rho*eta/Sigmav - 0.5;

    v1 = dt*1i*(c1 + v*c2);
    v2 = dt*.5*(v - Rho^2/denomQ*(v - vmin).*(vmax - v));
    fv = (1i*dxi*Rho/Sigmav)*v;
end

EXP_A = get_SV_matrix_exponential( Q, dt, xi, v1, v2, fv, psi_J, m_0, N );

%%%%////////////////////////////////////////////////////////
%%%% Construct Toepliz Array Of Arrays
%%%%////////////////////////////////////////////////////////
a2    = a^2;  
Cons2 = 24*a2*exp(-r*dt)/N;
zmin  = (1 - K)*dx;  %Kbar corresponds to zero

xi     = dxi*(1:N-1); %REDFINED FROM ABOVE
hvec = exp(-1i*zmin*xi).*(sin(xi/(2*a))./xi).^2./(2+cos(xi/a));

BETA = zeros(N,m_0,m_0);  %to access the (j,k)th toeplitz data, use BETA(:,j,k)
grand = zeros(1,N-1);
      
%%% NOTE the (k,j) rather than (j,k)
for j=1:m_0
    for k = 1:m_0
        for n=1:N-1  % MAKE this a .* if possible
            grand(n) = hvec(n)* EXP_A(k,j,n+1);  %n+1 since EXP_A has the first element defined
        end
        beta  = Cons2*real(fft([EXP_A(k,j,1)/(24*a2) grand]));   %%%%  NOTE: all toep matrices incorporate exp(-r*dt)
        toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];   
        BETA(:,j,k) = fft(toepM);
    end
end

%%%%////////////////////////////////////////////////////////
%%%% Initialize Continuation Value
%%%%////////////////////////////////////////////////////////
% CONT = repmat(p(1:K),1,m_0);  %continuation value in each state (initialize with all equal)
CONT = zeros(K,m_0);
ThetTemp = fft([ThetM(1:K);zeros(K,1)]);
for j=1:m_0
    for k = 1:m_0
        p = ifft(BETA(:,j,k).*ThetTemp);  
        CONT(:,j) = CONT(:,j)+ p(1:K);
    end
end

%%%%////////////////////////////////////////////////////////
%%%% LOOP through time
%%%%////////////////////////////////////////////////////////

kstr_vecInit = nbar*ones(m_0,1);   %in general one per state

for m = M-2:-1:0
    kstr_vec = kstr_vecInit;  %%%% FOR NOW!!! RESET EACH TIME (in future, we can start from previous known value to save cost)
    
    %Step 1: update THETA
    for j=1:m_0
        while kstr_vec(j)>2 && CONT(kstr_vec(j),j)> Gs(kstr_vec(j))
            kstr_vec(j) = kstr_vec(j) -1;
        end
        if kstr_vec(j)>=2
            xkstr = xmin +(kstr_vec(j) -1)*dx;

            Ck1 = CONT(kstr_vec(j)-1,j); Ck2 = CONT(kstr_vec(j),j); Ck3 = CONT(kstr_vec(j)+1,j);
            %%% Linear interp of payoff
            Gk2 = Gs(kstr_vec(j)); Gk3 = Gs(kstr_vec(j)+1);
            %%% xstr
            tmp1 = Ck2-Gk2; tmp2 = Ck3 - Gk3;
            xstrs  = ((xkstr+dx)*tmp1 - xkstr*tmp2)/(tmp1-tmp2);
        else
            kstr_vec(j) = 1; xstrs = xmin; xkstr = xmin; 
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
        
        THET(1:kstr_vec(j) -1,j) =  ThetM(1:kstr_vec(j)-1);
        Ck4 = CONT(kstr_vec(j) +2,j);
        
        THET(kstr_vec(j),j) = W*(.5 + dbar_0) -S_0*exp(xkstr)*(varthet_m10 +d_0)...
             + zeta4/8*(Ck1 -2*Ck2 + Ck3) +zeta3/3*(Ck2-Ck1)...
             + zeta2/4*(Ck1 +2*Ck2 - Ck3) - zeta*Ck2...
             -Ck1/24 +5/12*Ck2 +Ck3/8;
     
        THET(kstr_vec(j)+1,j) = W*dbar_1 - S_0*exp(xkstr+dx)*d_1    + zeta4/8*(-Ck2 +2*Ck3 - Ck4)...
                      + zeta3/6*(3*Ck2 - 4*Ck3 + Ck4) - .5*zeta2*Ck2...
                      + (Ck2 +10*Ck3 + Ck4)/12;
        THET(kstr_vec(j)+2:K-1,j) = (CONT(kstr_vec(j)+1:K-2,j)+10*CONT(kstr_vec(j)+2:K-1,j)+CONT(kstr_vec(j)+3:K,j))/12;
        THET(K,j)       = (13*CONT(K,j)+15*CONT(K-1,j)-5*CONT(K-2,j)+CONT(K-3,j))/48;

    end
    
    CONT = zeros(K,m_0);
    for k = 1:m_0
        ThetTemp = fft([THET(1:K,k);zeros(K,1)]);
        for j = 1:m_0
            p = ifft(BETA(:,j,k).*ThetTemp);
            CONT(:,j) = CONT(:,j)+ p(1:K);
        end
    end
    
end

%%%%////////////////////////////////////////////////////////
%%%% Interpolate to find price at v0
%%%%////////////////////////////////////////////////////////
k_0 = 2;  %k_0 will satisfy v(k_0) <= v0 < v(k_0 +1)
while v0 > v(k_0) && k_0 < m_0
    k_0 = k_0+1;
end
k_0 = k_0 - 1;

%%% Cubic spline
% k_int = [(k_0-2) (k_0-1) k_0 (k_0+1) (k_0+2)];
% v_int = [v(k_int(1)) v(k_int(2)) v(k_int(3)) v(k_int(4)) v(k_int(5))];
% Vals_int = [CONT(nnot,(k_int(1))) CONT(nnot,(k_int(2))) CONT(nnot,(k_int(3))) CONT(nnot,(k_int(4))) CONT(nnot,(k_int(5)))];
% Vals_int = max(Vals_int, Gs(nnot));
% 
% price = spline(v_int,Vals_int,v0);

%%% Linear Interpolation
%%% NOTE: we have assumed that are not in Case II (see paper.. this case occurs if 0 < ln(W/S_0) < Delta), otherwise we need 2 interpolations
Vals_Interp = [max(CONT(nnot,(k_0)),Gs(nnot)) max(CONT(nnot,(k_0+1)),Gs(nnot))];
price = Vals_Interp(1) + (Vals_Interp(2)-Vals_Interp(1))*(v0 - v(k_0))/(v(k_0+1)-v(k_0));
end

