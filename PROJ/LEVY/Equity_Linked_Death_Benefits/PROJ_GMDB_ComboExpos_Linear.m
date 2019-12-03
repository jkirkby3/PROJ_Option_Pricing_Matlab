function price = PROJ_GMDB_ComboExpos_Linear( P, Pbar, S_0, W, call, r, params_levy, params_mort, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Gauranteed Minimum Death Benefit Options using PROJ method
%        This version assumes a comination of exponentials model of mortality (see first reference below)
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
%
% Author: Zhimin Zhang  (Original Code)
%         Justin Lars Kirkby (Convert into common framework)
%
% References:  (1) Valuing Equity-Linked Death Benefits in General Exponential
%               Levy Models, J. Comput. and Appl. Math. 2019 (Z. Zhang, Y. Yong, W. Yu)
%              (2) Efficient Option Pricing By Frame Duality with The Fast
%               Fourier Transform, SIAM J. Financial Math., 2015 (J.L. Kirkby)
%
% ----------------------
% Contract/Model Params 
% ----------------------
% S_0 = initial stock price (e.g. 100)
% W   = strike  (e.g. 100)
% r   = interest rate (e.g. 0.05)
% T   = time remaining until maturity (in years, e.g. T=10)
%       NOTE: set T=-1 to price a perpetual contract (no expiry)
% call  = 1 for call (else put)
% params_levy = parameters of Levy Model
% params_mort = mortality params
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% P  = resolution parameter (increase P to use more basis elements)
% Pbar = gridwidth parameter (increase Pbar to increase the truncated density support)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Contract / Numerical Params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk=log(W/S_0);

a=2^P;
hat_a=2^Pbar;

N=hat_a*a;
Delta=1/a;Delta_s=2*pi/hat_a;

% Deterime if there is a Time to expiry in contract, else it's perpetual
if nargin < 8
    T = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse Death Model Params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = params_mort.lambda;  % expo rates
A = params_mort.A; % coefficients in expo combo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Levy Parse Model Params, set model inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = params_levy.model;

if model == 1 %BSM (Black Scholes Merton)
    sigma = params_levy.sigmaBSM;
    
    Psi=@(x)-1/2*sigma^2*x.^2;
    mu=r-Psi(-1i);
    
    phi=@(s)mu*1i*s-1/2*sigma^2*s.^2;
    
    if T > 0
        diffFgc=@(c)sum(A.*lambda*(mu+sigma^2*c)./(r+lambda-phi(-1i*c)).^2.*(1-exp(-(r+lambda-phi(-1i*c))*T)-T*exp(-(r+lambda-phi(-1i*c))*T).*(r+lambda-phi(-1i*c))));
    else
        diffFgc=@(c)sum(A.*lambda*(mu+sigma^2*c)./(r+lambda-phi(-1i*c)).^2);
    end
 
elseif model == 3 %NIG
    alpha = params_levy.alpha;
    beta = params_levy.beta;
    NIG_delta = params_levy.delta;
    sigma = params_levy.sigma;
    
    gamma=sqrt(alpha^2-beta^2);
    
    Psi=@(x)-1/2*sigma^2*x.^2-NIG_delta*(sqrt(alpha^2-(beta+1i*x).^2)-gamma);
    mu=r-Psi(-1i);
    
    phi=@(s)mu*1i*s-1/2*sigma^2*s.^2-NIG_delta*(sqrt(alpha^2-(beta+1i*s).^2)-gamma);
    
    if T > 0
        diffFgc=@(c)sum(A.*lambda*(mu+NIG_delta*sqrt(alpha^2-(beta+c).^2).*(beta+c))./(r+lambda-phi(-1i*c)).^2.*(1-exp(-(r+lambda-phi(-1i*c))*T)-T*exp(-(r+lambda-phi(-1i*c))*T).*(r+lambda-phi(-1i*c))));
    else
        diffFgc=@(c)sum(A.*lambda*(mu+NIG_delta*sqrt(alpha^2-(beta+c).^2).*(beta+c))./(r+lambda-phi(-1i*c)).^2);
    end
    
elseif model == 4 %MJD (Merton Jump Diffusion)
    sigma = params_levy.sigma;
    lambda_J = params_levy.lam;
    mu_J = params_levy.muj;
    sigma_J = params_levy.sigmaj;
    
    Psi=@(x)-sigma^2/2*x.^2+lambda_J*(exp(1i*x*mu_J-sigma_J^2/2*x.^2)-1);
    mu=r-Psi(-1i);

    phi=@(s)1i*s*mu-sigma^2/2*s.^2+lambda_J*(exp(1i*s*mu_J-sigma_J^2/2*s.^2)-1);
    
    if T > 0
        diffFgc=@(c)sum(A.*lambda*(mu+sigma^2*c+lambda_J*exp(mu_J*c+1/2*sigma_J^2*c).*(mu_J+sigma_J^2*c))./(r+lambda-phi(-1i*c)).^2.*(1-exp(-(r+lambda-phi(-1i*c))*T)-T*exp(-(r+lambda-phi(-1i*c))*T).*(r+lambda-phi(-1i*c))));
    else
        diffFgc=@(c)sum(A.*lambda*(mu+sigma^2*c+lambda_J*exp(mu_J*c+1/2*sigma_J^2*c).*(mu_J+sigma_J^2*c))./(r+lambda-phi(-1i*c)).^2);    
    end
    
elseif model == 5 %Kou Double Expo
    
    sigma = params_levy.sigma;
    lam_pois = params_levy.lam;   
    p = params_levy.p_up;    omega=p*lam_pois;  nv=(1-p)*lam_pois;
    v = params_levy.eta1;
    w = params_levy.eta2;
    
    D=1/2*sigma^2;

    %risk-neutral condition
    Psi=@(x)-D*x.^2-nv.*1i.*x./(v+1i.*x)+omega.*1i.*x./(w-1i.*x);
    mu=r-Psi(-1i);
    
    phi=@(s)mu.*1i.*s-D.*s.^2-nv.*1i.*s./(v+1i.*s)+omega.*1i.*s./(w-1i.*s);
    
    if T > 0
        diffFgc=@(c)sum(A.*lambda*(mu+sigma^2*c-nv*v./(v+c).^2+omega*w./(w-c).^2)./(r+lambda-phi(-1i*c)).^2.*(1-exp(-(r+lambda-phi(-1i*c))*T)-T*exp(-(r+lambda-phi(-1i*c))*T).*(r+lambda-phi(-1i*c))));
    else
        diffFgc=@(c)sum(A.*lambda*(mu+sigma^2*c-nv*v./(v+c).^2+omega*w./(w-c).^2)./(r+lambda-phi(-1i*c)).^2);
    end
    
elseif model == 8 % Variance Gamma
    sigma = params_levy.sigmaGBM; % geometric brownian motion add-on
    VG_mu = params_levy.theta;
    VG_sigma = params_levy.sigma;
    nv = params_levy.nu;
    
    Psi=@(x)-1/2*sigma^2*x.^2-1/nv*log(1-1i*nv*VG_mu*x+nv*VG_sigma^2/2*x.^2);
    mu=r-Psi(-1i);
    
    phi=@(s)mu*1i*s-1/2*sigma^2*s.^2-1/nv*log(1-1i*nv*VG_mu*s+nv*VG_sigma^2/2*s.^2);

    if T > 0
        diffFgc=@(c)sum(A.*lambda*(mu+sigma^2*c-(VG_mu*1i*c-VG_sigma^2*c)./(1-nv*VG_mu*c-1/2*nv*VG_sigma^2*c.^2))./(r+lambda-phi(-1i*c)).^2.*(1-exp(-(r+lambda-phi(-1i*c))*T)-T*exp(-(r+lambda-phi(-1i*c))*T).*(r+lambda-phi(-1i*c))));
    else
        diffFgc=@(c)sum(A.*lambda*(mu+sigma^2*c-(VG_mu*1i*c-VG_sigma^2*c)./(1-nv*VG_mu*c-1/2*nv*VG_sigma^2*c.^2))./(r+lambda-phi(-1i*c)).^2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


delta1=@(j,k)(j==k);

S=(0:(N-1))*Delta_s;
v=@(j)1-(delta1(j,1)+delta1(j,N))/2;
Fphi=@(s)12.*sin(s./2).^2./(s.^2.*(2+cos(s)));
Fphi1=Fphi(S/a); Fphi1(1)=1;

%coefs
if T > 0
    Ak=sum(repmat(A.*lambda,length(S),1)./(r+repmat(lambda,length(S),1)-repmat(phi(S).',1,2)).*(1-exp(-(r+repmat(lambda,length(S),1)-repmat(phi(S).',1,2))*T)),2);
    Ck=sum(repmat(A.*lambda,length(S),1)./(r+repmat(lambda,length(S),1)-repmat(phi(S-1i).',1,2)).*(1-exp(-(r+repmat(lambda,length(S),1)-repmat(phi(S-1i).',1,2))*T)),2);

    Fgc=@(c)sum(A.*lambda./(r+lambda-phi(-1i*c)).*(1-exp(-(r+lambda-phi(-1i*c))*T)));
else
    Ak=sum(repmat(A.*lambda,length(S),1)./(r+repmat(lambda,length(S),1)-repmat(phi(S).',1,2)),2);
    Ck=sum(repmat(A.*lambda,length(S),1)./(r+repmat(lambda,length(S),1)-repmat(phi(S-1i).',1,2)),2);

    Fgc=@(c)sum((A.*lambda)./(r+lambda-phi(-1i*c)));
end

% Compute Projections
x10=real(double(diffFgc(0))/Fgc(0))-Delta*(N/2-1);
X0=x10+(0:(N-1))*Delta;
beta0=a^(-1/2)/pi*real(fft((Ak.').*Fphi1.*v(1:N).*Delta_s.*exp(-1i*x10.*S)));

x11=real(double(diffFgc(1))/Fgc(1))-Delta*(N/2-1);
X1=x11+(0:(N-1))*Delta;
beta1=a^(-1/2)/pi*real(fft((Ck.').*Fphi1.*v(1:N).*Delta_s.*exp(-1i*x11.*S)));

% Compute Final Payoff and Value
nn=find((kk<X0)==1);
n=nn(1);
PHi0=zeros(1,length(S));

nn1=find((kk<X1)==1);
n1=nn1(1);
PHi1=zeros(1,length(S));

if call == 1
    PHi0(n)=a^(1/2)*(X0(n)-kk)-a^(3/2)/2*(kk-X0(n))^2+1/2*a^(-1/2);
    PHi0(n-1)=a^(1/2)*(X0(n-1)+1/a-kk)-a^(3/2)/2*(1/(a^2)-(kk-X0(n-1))^2);
    PHi0(n+1:length(S))=1/sqrt(a);
    
    PHi1(n1)=a^(1/2)*(X1(n1)-kk)-a^(3/2)/2*(kk-X1(n1))^2+1/2*a^(-1/2);
    PHi1(n1-1)=a^(1/2)*(X1(n1-1)+1/a-kk)-a^(3/2)/2*(1/(a^2)-(kk-X1(n1-1))^2);
    PHi1(n1+1:length(S))=1/sqrt(a);
    
    price=S_0*beta1*PHi1.'-W*beta0*PHi0.';
else
    PHi0(n)=a^(1/2)*(kk-X0(n)+1/a)-a^(3/2)/2*(1/(a^2)-(kk-X0(n))^2);
    PHi0(n-1)=a^(1/2)*(kk-X0(n-1))-a^(3/2)/2*(kk-X0(n-1))^2+1/2*a^(-1/2);
    PHi0(1:n-2)=1/sqrt(a);

    PHi1(n1)=a^(1/2)*(kk-X1(n1)+1/a)-a^(3/2)/2*(1/(a^2)-(kk-X1(n1))^2);
    PHi1(n1-1)=a^(1/2)*(kk-X1(n1-1))-a^(3/2)/2*(kk-X1(n1-1))^2+1/2*a^(-1/2);
    PHi1(1:n1-2)=1/sqrt(a);
    
    price=W*beta0*PHi0.'-S_0*beta1*PHi1.';
end

end

