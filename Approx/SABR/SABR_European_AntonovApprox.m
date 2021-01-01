function price = SABR_European_AntonovApprox(f,K,T,call,r,ModParams,UB)
% About: Price European Option under SABR model using Antonov approximation
% K: strike
% f: initial forward value
% call: 1 if call option, else put-call
% ModParams: container of model parameters (.beta, .alpha, .v0, .rho)
% UB: upper bound numerical parameter (optional)

if nargin < 7
   UB = 100; 
end

beta = ModParams.beta;
alpha = ModParams.alpha;
v0 = ModParams.v0;
rho = ModParams.rho;

mimic_parameters = mimicprocess(f,v0,beta,rho,alpha,K,T);                              % Mimic parameters
price = Vcallapprox(f,mimic_parameters(1),mimic_parameters(2),mimic_parameters(3),K,T,UB); % Call price;
price = exp(-r*T)*price;

if call ~= 1   % use put-call parity to price
    price = price - exp(-r*T)*(f - K);
end
end

function Vcall = Vcallapprox(f,v0,beta,nu,K,T,UB)
% Function to calculate the forward value of an European call with
% Antonov's method in the zero correlation case and with G approximated by
% (11)

%f     = initial forward value
%v0    = initial volatility 
%beta  = exponent
%nu    = vol-vol
%correlation rho is assumed to be zero

q     = K^(1-beta)/(1-beta);
q0    = f^(1-beta)/(1-beta);
smin  = asinh(nu*abs(q-q0)/v0);
splus = asinh(nu*(q+q0)/v0);
eta   = abs(1/(2*(beta-1)));

g  = @(s)s.*coth(s)-1;
R  = @(t,s)1+3*t*g(s)./(8*s.^2)-5*t^2*(-8*s.^2+3*g(s).^2+24*g(s))./(128*s.^4)+35*t^3*(-40*s.^2+3*g(s).^3+24*g(s).^2+120*g(s))./(1024*s.^6);
dR = @(t,s)exp(t/8)-(3072+384*t+24*t^2+t^3)/3072;
G  = @(t,s)sqrt(sinh(s)./s).*exp(-s.^2/(2*t)-t/8).*(R(t,s)+dR(t,s));

phi = @(s)2*atan(sqrt((sinh(s).^2-sinh(smin).^2)./(sinh(splus).^2-sinh(s).^2)));
psi = @(s)2*atanh(sqrt((sinh(s).^2-sinh(splus).^2)./(sinh(s).^2-sinh(smin).^2)));

Vcall = max(0,f-K)+2/pi*sqrt(K*f).*(quad(@(s)sin(eta*phi(s))./(sinh(s)).*G(T*nu^2,s),smin,splus)+sin(eta*pi)*quad(@(s)exp(-eta*psi(s))./sinh(s).*G(T*nu^2,s),splus,UB));
end


function parameters = mimicprocess(f,v0,beta,rho,nu,K,T)
% Function to find the mimic parameters to find the forward value of an European call with
% Antonov's method in the correlated case

%f     = initial forward value
%v0    = initial volatility 
%beta  = exponent
%nu    = vol-vol
%K     = Strike price
%T     = Time to maturity
%correlation rho is assumed to be non-zero

%alpha_mimic = initial volatility of the mimic process
%beta_mimic = exponent of the mimic process
%nu_mimic = vol-vol of the mimic process

beta_mimic  = beta;                                                             % (21)
nu_mimic    = sqrt(nu^2-3*(nu^2*rho^2+v0*nu*rho*(1-beta)*f^(beta-1))/2);     % (22)
dq          = (K^(1-beta)-f^(1-beta))/(1-beta);                                 % Above (14)
dq_mimic    = (K^(1-beta_mimic)-f^(1-beta_mimic))/(1-beta_mimic);          % Above (14)
vmin        = sqrt(nu^2*dq^2+2*rho*nu*dq*v0+v0^2);                            % Above (14)
Phi         = ((vmin+rho*v0+nu*dq)/((1+rho)*v0))^(nu_mimic/nu);           % Between (19) and (20)
v_mimic0    = 2*Phi*dq_mimic*nu_mimic/(Phi^2-1);                               % (19)
u0          = (dq*nu*rho+v0-vmin)/(dq*nu*sqrt(1-rho^2));                     % Between (14) and (15)
L           = vmin*(1-beta)/(K^(1-beta)*nu*sqrt(1-rho^2));                      % Between (14) and (15)

if L<1                                                                          % Between (14) and (15)
    I    = 2*(atan((u0+L)/sqrt(1-L^2))-atan(L/sqrt(1-L^2)))/sqrt(1-L^2);
else
    I    = log((u0*(L+sqrt(L^2-1))+1)/(u0*(L-sqrt(L^2-1))+1))/sqrt(L^2-1);
end 

phi0            = acos(-(dq*nu+v0*rho)/vmin);                                % Between (14) and (15)
betamin         = -beta*rho*(pi-phi0-acos(rho)-I)/(2*(1-beta)*sqrt(1-rho^2));   % Between (14) and (15)
v_mimic1         = v_mimic0*nu_mimic^2*(((beta-beta_mimic)*log(K*f)+log(v0*vmin)-log(v_mimic0*sqrt(dq_mimic^2*nu_mimic^2+v_mimic0^2)))/2-betamin)/((Phi^2-1)*log(Phi)/(Phi^2+1)); %(20)
alpha_mimic     = v_mimic0+T*v_mimic1; %(18)

parameters = [alpha_mimic,beta_mimic,nu_mimic];
end