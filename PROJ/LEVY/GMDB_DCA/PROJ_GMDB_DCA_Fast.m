function [price, opt, L] = PROJ_GMDB_DCA_Fast(proj_params, S_0, gmdb_params, r, q, modelInput)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for DCA-Style Garuanteed Minimum Withdraw Benefit (GMWB) using PROJ method
%        This version is based on a dollar cost average style investment account (see reference below)
%
% Terminal Payoff:  Payoff(tau) = L*exp(g*tau) + (Gam(tau) - L*exp(g*tau))^+
%                      Gam(tau) = S_M * sum_{m=0}^M(alpha*gamma / S_m)
%                          tau  = time of death (discrete periods)
%                            M  = number of periods until time of death (each period length dt)
%
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
%
% Author: Justin Lars Kirkby
% References: 1) Equity-Linked  Guaranteed Minimum Death Benefits with Dollar Cost Averaging, J.L.Kirkby & D.Nguyen, 2021
%
% ----------------------
% Contract/Model Params 
% ----------------------
% S_0 = initial stock price (e.g. 100)
% r   = interest rate (e.g. 0.05)
% q   = dividend yield (e.g. 0.05)
% M   = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
% gmdb_params = container of GMDB contract params, see below
% modelInput =  model inputs, see below
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% proj_params = numerical params
%   proj_params.N = number of basis elements, e.g. N = 2^10
%   proj_params.L1 = gridwidth param, e.g. L1 = 8
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------
% GMDB Contract Params
% ------------------
L = gmdb_params.L;   % Guarantee Level: Set L = -1 to use ATMF value for L
alpha = gmdb_params.alpha;  % Period premium payment, paid every dt time units
gamma = gmdb_params.gamma;  % Proportion of investment retained by policyholder (fee is 1-gamma)
contract_type = gmdb_params.contract_type;  % Contract type: 1 = GMDB, 2 = GMDB-RS (Ratchet strike)
p = gmdb_params.death_prob;  % death probability distribution, must be consistent with dt
g = gmdb_params.g;

% ------------------
% Model Inputs
% ------------------
dt = modelInput.dt;  % Time increment, premiums paid / underlying S is monitored every dt
phiR = modelInput.rnCHF;  % Risk neutral CHF for time period dt


call = 1;
ER = 0;
  
Z = gen_func(-r, dt, p);
if g == 0
    Zrg = Z;
else
    Zrg = gen_func(-(r-g), dt, p);
end

if L == -1
    MF = gen_func(r - q - g, dt, p);
    Zg = gen_func(-g, dt, p);
    L = alpha * gamma * (exp((r-q)*dt)*MF - Zg) / (exp((r-q)*dt) - 1);
end

Mmax = length(p); Tmax = Mmax*dt;
pr_alpha = getTruncationAlpha(Tmax, proj_params.L1, modelInput, modelInput.model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = proj_params.N;
dx = 2*pr_alpha/(N-1); a = 1/dx;
A    = 32*a^4;
C_aN = A/N;



%%% SHIFTS
x1    = zeros(1,Mmax);

if contract_type == 1
    strikes = S_0*L*exp(g*dt*(1:Mmax)) ./ (alpha*gamma*(2:Mmax + 1));
else
    strikes = S_0*ones(1,Mmax);
end


for m=1:Mmax
    if m == 1
        x1(m) = ER;
    else
        x1(m) = ER + log(1+exp(x1(m-1)));  %%BENHAMOU SHIFT
    end
    %x1(m) = log(m) + .5*(m+1)*ER;    %% LOWER BOUND SHIFT derived in APROJ paper
end

Nm   = floor(a*(x1-ER));
x1   = ER + (1-N/2)*dx + Nm*dx;
NNM  = N + Nm(Mmax-1);   %Number of columns of PSI

%%% Now check that we wont fall off the grid later
for m=1:Mmax
    ystar = log((m+1)*strikes(m)/S_0 -1);
    nbar  = floor((ystar-x1(m))*a+1);
    % fprintf('%.0f\n', nbar);
    if nbar + 1 > N  
        proj_params.L1 = proj_params.L1*1.25;
        [price, opt, L] = PROJ_GMDB_DCA_Fast(proj_params, S_0, gmdb_params, r, q, modelInput);
        return;
    end
end

dxi   = 2*pi*a/N;
xi    = dxi*(1:(N-1))';
PhiR = [1; phiR(xi)];

b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
zeta  = (sin(xi/(2*a))./xi).^4./(b0 + b1*cos(xi/a) +b2*cos(2*xi/a) +b3*cos(3*xi/a));
AA = 1/A;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PSI Matrix: 5-Point GAUSSIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI = make_PSI(N,NNM,x1(1),dx,dxi);

%%%%%%%%%%%%%%
% STEP 1) Value the European!!!
%%%%%%%%%%%%%%
beta  = [AA; zeta.*PhiR(2:N).*exp(-1i*x1(1)*xi)];   %grand(end)=.5*grand(end);
beta  = real(fft(beta));

s = european_price(zeta,PhiR,xi,dx, r, q, dt, strikes(1), S_0, ER, N, a, call);
s = p(1) * 2 * s;

%%%%%%%%%%%%%%
% STEP 2) Value the rest
%%%%%%%%%%%%%%

PhiR  = C_aN*PhiR;
beta  = PSI(:,1:N)*beta.*PhiR;  %Nm(1)=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Loop to find PSI_M
for n=2:Mmax
    opt_v = intermediate_asian_price(N, dx, dt, xi, zeta, beta, x1(n), n, r, q, strikes(n), S_0);
    % fprintf('%.12f \n',opt_v);
    
    beta(2:N) = zeta.*beta(2:N).*exp(-1i*x1(n)*xi); beta(1) = AA;
    beta      = real(fft(beta));
    if n < Mmax
        beta      = PSI(:,Nm(n)+1:Nm(n)+N)*beta.*PhiR;
    end
    
    s = s + p(n) * (n + 1) * opt_v; 
end


opt = s * alpha * gamma / S_0;
price = L*Zrg - alpha * (exp(r*dt) - Z) / (exp(r*dt) - 1) + opt;

end


function PSI = make_PSI(N,NNM,x_1,dx,dxi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PSI Matrix: 5-Point GAUSSIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI     = zeros(N,NNM);    %The first row will remain ones
PSI(1,:) = ones(1,NNM);

%%%% Sample
Neta  = 5*(NNM) + 15;   %sample size
Neta5 = (NNM) + 3;
g2    = sqrt(5-2*sqrt(10/7))/6;
g3    = sqrt(5+2*sqrt(10/7))/6;
v1    = .5*128/225; 
v2    = .5*(322+13*sqrt(70))/900; 
v3    = .5*(322 - 13*sqrt(70))/900;


thet                 = zeros(1,Neta);   %sample initialized
thet(5*(1:Neta5)-2)  = x_1 -1.5*dx + dx*(0:Neta5-1);
thet(5*(1:Neta5)-4)  = x_1 -1.5*dx + dx*(0:Neta5-1) - dx*g3;
thet(5*(1:Neta5)-3)  = x_1 -1.5*dx + dx*(0:Neta5-1) - dx*g2;
thet(5*(1:Neta5)-1)  = x_1 -1.5*dx + dx*(0:Neta5-1) + dx*g2;
thet(5*(1:Neta5))    = x_1 -1.5*dx + dx*(0:Neta5-1) + dx*g3;


%%%% Weights
sig      = [-1.5-g3, -1.5-g2, -1.5, -1.5+g2, -1.5+g3, -.5-g3, -.5-g2, -.5, -.5+g2, -.5+g3,];
sig(1:5) = (sig(1:5) + 2).^3/6;
sig(6:10) = 2/3 - .5*(sig(6:10)).^3 - (sig(6:10)).^2;

sig([1 5 6 10]) = v3*sig([1 5 6 10]); 
sig([2 4 7 9]) = v2*sig([2 4 7 9]); 
sig([3 8]) = v1*sig([3 8]);

%%%% Fill Matrix
zz  = exp(1i*dxi*log(1+exp(thet)));
thet   = zz; 


for j=2:N-1
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
end

function Val = intermediate_asian_price(N, dx, dt, xi, zeta, chf, x_1, M, r, q, W, S_0)
call = 1;
a = 1/dx;
A    = 32*a^4;
AA = 1/A;
C_aN = A/N;
T = M*dt;

%%%%% FINAL VALUE
ystar = log((M+1)*W/S_0 -1);
nbar  = floor((ystar-x_1)*a+1);
C     = S_0/(M+1);
D     = W - C;
x_1   = ystar- (nbar-1)*dx;

beta(2:N) = zeta.*chf(2:N).*exp(-1i*x_1*xi); beta(1)=AA;
beta      = real(fft(beta)).';



Cc1 = C*( exp(-7/4*dx)/54 + exp(-1.5*dx)/18 + exp(-1.25*dx)/2 + 7*exp(-dx)/27 )/20;

Cc2 = C*.05*(28/27 + exp(-7/4*dx)/54 + exp(-1.5*dx)/18 + exp(-1.25*dx)/2 + 14*exp(-dx)/27 ...
              + 121/54*exp(-.75*dx) + 23/18*exp(-.5*dx) + 235/54*exp(-.25*dx));
          
Cc3 = C*( (28 + 7*exp(-dx))/3 ...
              + ( 14*exp(dx) + exp(-7/4*dx) + 242*cosh(.75*dx) + 470*cosh(.25*dx))/12 ...
              +.25*(exp(-1.5*dx) + 9*exp(-1.25*dx) + 46*cosh(.5*dx)))/90;
          
Cc4 = C*( 14/3*(2+cosh(dx)) ...
              + .5*(cosh(1.5*dx) + 9*cosh(1.25*dx) +23*cosh(.5*dx))...
              +  1/6*(cosh(7/4*dx) + 121*cosh(.75*dx) +235*cosh(.25*dx)))/90;

G           = zeros(nbar+1,1);
E           = exp(ystar-(nbar-1)*dx+dx*(0:nbar));

G(nbar+1)   = D/24    - Cc1*E(nbar+1);
G(nbar)     = .5*D    - Cc2*E(nbar);
G(nbar-1)   = 23*D/24 - Cc3*E(nbar-1);
G(1:nbar-2) = D       - Cc4*E(1:nbar-2); 

Val = C_aN*exp(-r*T)*sum(beta(1:nbar+1).*G);
if call==1  %Call Option, use Put-Call-Parity
    if r - q == 0
        mult = M + 1;
    else
        mult = (exp((r-q)*T*(1+1/M))-1)/(exp((r-q)*dt)-1);
    end
    Val = Val + C*exp(-r*T)*mult - W*exp(-r*T);
end
Val = max(0, Val);

end


function price = european_price(zeta,PhiR,xi,dx, r, q, T, W, S_0, c1, N, a, call)

W = 2*W - S_0;  % we shift the strike then divide by 2 at the end, this is asian with M=2

lws = log(W/S_0);
lam = c1 -(N/2 -1)*dx;
nbar = floor(a*(lws-lam)+1);
if nbar>=N
    nbar = N-1;
end
xmin = lws - (nbar-1)*dx;


Cons = 32*a^4;

beta  = [1/Cons; zeta.*PhiR(2:N).*exp(-1i*xmin*xi)];   %grand(end)=.5*grand(end);
beta  = real(fft(beta)).';

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
            
if call == 1  % Use put-call parity
    price = Cons*exp(-r*T)/N*G*(beta(1:length(G))') + S_0*exp(-q*T) - W*exp(-r*T);
else
    price = Cons*exp(-r*T)/N*G*(beta(1:length(G))');
end

price = 0.5*max(price, 0);  % Protect against deep out of money case

end
