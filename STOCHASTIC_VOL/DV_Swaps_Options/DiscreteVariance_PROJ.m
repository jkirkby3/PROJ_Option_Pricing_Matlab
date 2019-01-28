function price = DiscreteVariance_PROJ( N,alph,M,r,T,K,m_0,psi_J,model, modparam, gridMethod, gamma, varGridMult, contract )
% N = #basis points
% alph = log-asset grid width param
% M = # Monitoring dates (not including S_0)
% r = interest rate
% T = time to maturity
% K = strike (only matters for an option, but is always required)
% m_0 = number of Variance States in markov chain approx
% psi_J = characteristic exponent of the jump component of model, function handle
% model:
%        1 = HESTON:      Sigmav, v0, rho, eta, theta
%        2 = STEIN-STEIN: Sigmav, v0, rho, eta, theta
%        3 = 3/2 MODEL:   Sigmav, v0, rho, eta, theta
%        4 = 4/2 MODEL:   Sigmav, v0, rho, eta, theta, aa, bb
%        5 = HULL-WHITE:  Sigmav, v0, rho
%        6 = SCOTT:       Sigmav, v0, rho, eta, theta
%        7 = ALPHA-HYPER: Sigmav, v0, rho, eta, theta
%
% modparam : parameters of the model, see above
% gridMethod - Always use with 4 for now (this is the nonuniform grid)
% gamma = variance gridwidth multiplier
% contract: 
%           1 = Variance Swap, 
%           2 = Volatility Swap, 
%           3 = Call on Variance, 
%           4 = Put on Variance
% NOTE: right now only 1 and 3 are supported
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx   = 2*alph/(N-1);
dt   = T/M;
a    = 1/dx;
A    = 32*a^4;
C_aN = A/N;
xmin = (1-N/2)*dx; 

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

    Rho = modparam.rho; Sigmav = modparam.Sigmav; v0 = modparam.v0; eta = modparam.eta; av = modparam.av; theta = modparam.theta;
    mu_func = @(v)eta - theta*exp(av*v);
    sig_func = @(v)Sigmav*(v>-100);
    
    %%% Estimate Mean and Variance Using Stein Stein and First order approx
    EtaBar = theta*av*exp(av*v0); ThetaBar = (eta - theta*exp(av*v0)*(1-av*v0))/EtaBar;
    mu_H = exp(-EtaBar*t)*v0 + ThetaBar*(1-exp(-EtaBar*t));
    sig2_H = Sigmav^2/(2*EtaBar)*(1 -exp(-2*EtaBar*t));
end

if model == 2 || model == 4  %For Stein-Stein and 4/2 more sensitive
    lx = max(0.01, mu_H - gamma*sqrt(sig2_H));  %variance grid lower bound
elseif model == 6 || model == 7   %For Scott/alpha-Hyper, allow lx negative
    lx = mu_H - gamma*sqrt(sig2_H);
else
    lx = max(0.00001,mu_H - gamma*sqrt(sig2_H));  %variance grid lower bound
end
ux = mu_H + gamma*sqrt(sig2_H);  %variance grid upper bound

tic
center = v0; %this is where grid clusters... we can experiment with other choices.. 
%%% [Q,v]  = General_Q_Matrix_DiscreteVar(m_0,mu_func,sig_func,lx,ux,gridMethod,center,varGridMult);
boundaryMethod = 1;
[Q,v]  = Q_Matrix_AllForms(m_0,mu_func,sig_func,lx,ux,gridMethod, varGridMult, center, boundaryMethod);

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
    c1 = Rho*theta/Sigmav;
    c2 = Rho*(eta/Sigmav + Sigmav/2); c3 = .5; c4 = r - psi_J(-1i);
    expv = exp(v); expv2 = expv.^2;
    v1 = dt*1i*(c1*exp((1+av)*v) - c2*expv - c3*expv2 + c4);
    v2 = dt*.5*(1-Rho^2)*expv2;
    fv = 1i*dxi*Rho/Sigmav*expv;
end


% Compute Matrix Exponentials for each xi(j)
EXP_A = get_SV_matrix_exponential( Q, dt, xi, v1, v2, fv, psi_J, m_0, N );


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% PSI Matrix: 5-Point GAUSSIAN
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
NNM = N;  %% in this case, same dimension
PSI = zeros(N-1,NNM);    %The first row (left undefined) will remain ones (since the chf will always be one at that point)

%%%% Sample
Neta  = 5*(NNM) + 15;   %sample size
Neta5 = (NNM) + 3;
g2    = sqrt(5-2*sqrt(10/7))/6;
g3    = sqrt(5+2*sqrt(10/7))/6;
v1    = .5*128/225; v2 = .5*(322+13*sqrt(70))/900;  v3 = .5*(322 - 13*sqrt(70))/900;

thet                 = zeros(1,Neta);   %sample initialized
thet(5*(1:Neta5)-2)  = xmin -1.5*dx + dx*(0:Neta5-1);
thet(5*(1:Neta5)-4)  = xmin -1.5*dx + dx*(0:Neta5-1) - dx*g3;
thet(5*(1:Neta5)-3)  = xmin -1.5*dx + dx*(0:Neta5-1) - dx*g2;
thet(5*(1:Neta5)-1)  = xmin -1.5*dx + dx*(0:Neta5-1) + dx*g2;
thet(5*(1:Neta5))    = xmin -1.5*dx + dx*(0:Neta5-1) + dx*g3;

%%%% Weights
sig      = [-1.5-g3, -1.5-g2, -1.5, -1.5+g2, -1.5+g3, -.5-g3, -.5-g2, -.5, -.5+g2, -.5+g3,];
sig(1:5) = (sig(1:5) + 2).^3/6;
sig(6:10) = 2/3 - .5*(sig(6:10)).^3 - (sig(6:10)).^2;

sig([1 5 6 10]) = v3*sig([1 5 6 10]); sig([2 4 7 9]) = v2*sig([2 4 7 9]); sig([3 8]) = v1*sig([3 8]);
sig = C_aN*sig;  

%%%% Fill Matrix  (NOTE: this can be made MORE EFFICIENT by using symmetery of x^2)
zz  = exp(1i*dxi*thet.^2); %% in general, 1i*dxh(thet)
thet   = zz; 

for j=1:N-1  %Note: first row is not ones anymore
    PSI(j,:) =  sig(1)*(thet(1:5:Neta-19) + thet(20:5:Neta)) ...
              + sig(2)*(thet(2:5:Neta-18) + thet(19:5:Neta-1)) ...
              + sig(3)*(thet(3:5:Neta-17)  + thet(18:5:Neta-2)) ...
              + sig(4)*(thet(4:5:Neta-16)  + thet( 17:5:Neta-3)) ...
              + sig(5)*(thet(5:5:Neta-15)  + thet( 16:5:Neta-4)) ...
              + sig(6)*(thet(6:5:Neta-14)  + thet( 15:5:Neta-5)) ...
              + sig(7)*(thet(7:5:Neta-13)  + thet( 14:5:Neta-6)) ...
              + sig(8)*(thet(8:5:Neta-12)  + thet( 13:5:Neta-7)) ...
              + sig(9)*(thet(9:5:Neta-11)  + thet( 12:5:Neta-8)) ...
              + sig(10)*(thet(10:5:Neta-10)  + thet( 11:5:Neta-9));

    thet = thet.*zz;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find phi_{Y_1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

xi    = dxi*(1:N-1)';  %REDEFINED from above, excludes the zero at xi_1
b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
zeta  = (sin(xi/(2*a))./xi).^4./(b0 + b1*cos(xi/a) +b2*cos(2*xi/a) +b3*cos(3*xi/a));
hvec  = exp(-1i*xmin*xi).*zeta;  %dont define first element

%Note: PHIY is not defined for xi = 0
PHIY_old = zeros(N-1,m_0);
PHIY_new = zeros(N-1,m_0);  %This serves intermediate storage for H,Beta,Phi_z

%Find beta in first stage (using chf of Phi_Y1)
for j = 1:m_0
   %Step 1: characteristic function of log return
   for n = 1:N-1
    PHIY_old(n,j) = sum(EXP_A(1:m_0,j,n+1));  %n+1 since EXP_A is defined for xi = 0
   end
   
   %Step 2: invert characteristic function of log return (ie this is beta)
   BetaTemp =  real(fft([1/A; PHIY_old(:,j).*hvec] )); %note: dont need PHIY_old(1,j)/A, since in this case PHIY_old(1,j) = 1
   
   %Step 3: Phi_{Y_1}^j 
   PHIY_new(:,j) = PSI*BetaTemp; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find PHI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
PHI  = ones(m_0,m_0,N-1); 
%%% use PHIY_old for temp storage
grand = zeros(N-1,1); %%%temp vec
for j = 1:m_0
    for k=1:m_0
        %First Invert chf to get p_{j,k}
        for n=1:N-1  % MAKE this a .*
            grand(n) = hvec(n)* EXP_A(k,j,n+1);  %n+1 since EXP_A has the first element defined
        end
        BetaTemp  = real(fft([EXP_A(k,j,1)/A; grand])); 
        PHI(j,k,:) = PSI*BetaTemp;
        %PHI(j,k,:) = PSI*real(fft([EXP_A(k,j,1)/A; grand]));
    end
end
clear EXP_A

for m =2:M
    for n = 1:N-1
        PHIY_new(n,:) = PHIY_new(n,:)*PHI(:,:,n).';
    end
end

%%%%////////////////////////////////////////////////////////
%%%% Interpolate to find bracketing initial volatilities
%%%%////////////////////////////////////////////////////////
k_0 = 2;  %k_0 will satisfy v(k_0) <= v0 < v(k_0 +1)
while v0 >= v(k_0) && k_0 < m_0
    k_0 = k_0+1;
end
k_0 = k_0 - 1;

%%%%%%%%%%%%%%%%%%%%%%%%% 
cubicTerminal = 1;  %Set to 1 to use cubic projection in final stage, else does linear 
%%%%%%%%%%%%%%%%%%%%%%%%% 
if contract == 1 %Variance Swap
    if cubicTerminal == 1
        xmin = -dx;
    else
        xmin = 0;
    end
elseif contract == 3 %Call on Variance
    if cubicTerminal == 1
        xmin = K*T - dx;
    else
        xmin = K*T;
    end
end

if cubicTerminal == 1
    if contract == 1 || contract == 3 %Variance Swap or Call on Variance   
        grid = -dx + dx*(0:N-1);  %NOTE: this holds for call on variance too, since we subtract out the KT

        grid(1) = grid(1)/24 + dx/20;
        grid(2) = dx*7/30;  %note: grid(2) = xmin + dx = 0;
        grid(3) = grid(3)*23/24 + dx/20;      
    end    
else  %Use Linear Projection at the end
    if contract == 1 || contract == 3 %Variance Swap or Call on Variance   
        grid = dx*(0:N-1);  %NOTE: this holds for call on variance too, since we subtract out the KT
        grid(1) = dx/6;
    end
    A = 24*a^2;
    C_aN  = A/N; 
    zeta = (sin(xi/(2*a))./xi).^2./(2+cos(xi/a));
end

if xmin ~= 0
    hvec = exp(-1i*xmin*xi).*zeta;
else
    hvec = zeta;
end

    
vals = [0 0];
ks = [k_0 k_0+1];    

for l = 1:2
    j = ks(l);
    BetaTemp  = real(fft([1/A; hvec.*PHIY_new(:,j)])); 
    if contract == 1 || contract == 3 %Variance Swap or call
        vals(l) = grid(1:N/2)*BetaTemp(1:N/2);
        vals(l) = C_aN*vals(l);
    end
end
toc

if gridMethod == 5 || gridMethod == 6  %Then NO Interpolation Needed
    Approx = vals(1);
else
    Approx = vals(1) + (vals(2)-vals(1))*(v0 - v(k_0))/(v(k_0+1)-v(k_0));
end
    
if contract == 1 %Variance Swap
    price = Approx/T;  %can generalize to include anualization factor if needed
elseif contract == 3 %Variance Call
    price = exp(-r*T)*Approx/T;
else
    fprintf('Only contract types 1 and 3 are currently supported \n')
end

end
