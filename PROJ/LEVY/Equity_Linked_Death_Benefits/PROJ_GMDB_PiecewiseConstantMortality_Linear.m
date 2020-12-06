function price = PROJ_GMDB_PiecewiseConstantMortality_Linear( P, Pbar, S_0, W, call, r, params_levy, params_mort)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Gauranteed Minimum Death Benefit Options using PROJ method
%        This version assumes a piecewise constant model of mortality (see first reference below)
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
% call  = 1 for call (else put)  ... NOTE: currently only a PUT contract is supported
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse Death Model Params (Piceweise constant forces of mortality)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qx = params_mort.qx;  % life table
x = params_mort.x; % current age
max_age = params_mort.max_age;

n_years=max_age-x; %n year contigent
px=1-qx;
mux=-log(px(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Levy Parse Model Params, set model inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = params_levy.model;

if model == 1 %BSM (Black Scholes Merton)
    sigma = params_levy.sigmaBSM;
    
    Psi=@(x)-1/2*sigma^2*x.^2;
    mu=r-Psi(-1i);
    
    phi=@(s)mu*1i*s-1/2*sigma^2*s.^2;
    

elseif model == 3 %NIG
    alpha = params_levy.alpha;
    beta = params_levy.beta;
    NIG_delta = params_levy.delta;
    sigma = params_levy.sigma;
    
    gamma=sqrt(alpha^2-beta^2);
    
    Psi=@(x)-1/2*sigma^2*x.^2-NIG_delta*(sqrt(alpha^2-(beta+1i*x).^2)-gamma);
    mu=r-Psi(-1i);
    
    phi=@(s)mu*1i*s-1/2*sigma^2*s.^2-NIG_delta*(sqrt(alpha^2-(beta+1i*s).^2)-gamma);
    
    
elseif model == 4 %MJD (Merton Jump Diffusion)
    sigma = params_levy.sigma;
    lambda_J = params_levy.lam;
    mu_J = params_levy.muj;
    sigma_J = params_levy.sigmaj;
    
    Psi=@(x)-sigma^2/2*x.^2+lambda_J*(exp(1i*x*mu_J-sigma_J^2/2*x.^2)-1);
    mu=r-Psi(-1i);

    phi=@(s)1i*s*mu-sigma^2/2*s.^2+lambda_J*(exp(1i*s*mu_J-sigma_J^2/2*s.^2)-1);
    
    
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
    
    
elseif model == 8 % Variance Gamma
    sigma = params_levy.sigmaGBM; % geometric brownian motion add-on
    VG_mu = params_levy.theta;
    VG_sigma = params_levy.sigma;
    nv = params_levy.nu;
    
    Psi=@(x)-1/2*sigma^2*x.^2-1/nv*log(1-1i*nv*VG_mu*x+nv*VG_sigma^2/2*x.^2);
    mu=r-Psi(-1i);
    
    phi=@(s)mu*1i*s-1/2*sigma^2*s.^2-1/nv*log(1-1i*nv*VG_mu*s+nv*VG_sigma^2/2*s.^2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1=mu-Delta*(N/2-1);
X=x1+(0:(N-1))*Delta;
delta1=@(j,k)(j==k);

S=(0:(N-1))*Delta_s;
v=@(j)1-(delta1(j,1)+delta1(j,N))/2;
Fphi=@(s)12.*sin(s./2).^2./(s.^2.*(2+cos(s)));
Fphi1=Fphi(S/a); Fphi1(1)=1;

%coefs
Ak=zeros(1, N);
for k=0:N-1
    Ak(k+1)=0;
    for j=0:n_years-1
        Ak(k+1)=Ak(k+1)+prod(px(x+1:x+j))*exp(j*mux(x+j+1))*mux(x+j+1)*(exp(-(r+mux(x+j+1)-phi(S(k+1)))*j)-exp(-(r+mux(x+j+1)-phi(S(k+1)))*(j+1)))/(r+mux(x+j+1)-phi(S(k+1)));
    end
end

Bk=zeros(1, N);
for k=0:N-1
    Bk(k+1)=0;
    for j=0:n_years-1
        Bk(k+1)=Bk(k+1)+prod(px(x+1:x+j))*exp(j*mux(x+j+1))*mux(x+j+1)*(exp(-(r+mux(x+j+1)-phi(S(k+1)-1i))*j)-exp(-(r+mux(x+j+1)-phi(S(k+1)-1i))*(j+1)))/(r+mux(x+j+1)-phi(S(k+1)-1i));
    end
end

beta0=a^(-1/2)/pi*real(fft((Ak).*Fphi1.*v(1:N).*Delta_s.*exp(-1i*x1.*S)));
beta1=a^(-1/2)/pi*real(fft((Bk).*Fphi1.*v(1:N).*Delta_s.*exp(-1i*x1.*S)));

nn=find((kk<X)==1);
n=nn(1);
PHi=zeros(1,length(S));

if call == 1
    fprintf('NOT IMPLEMENTED for call option yet');
    price = -123456789;
    
else
    PHi(n)=a^(1/2)*(kk-X(n)+1/a)-a^(3/2)/2*(1/(a^2)-(kk-X(n))^2);
    PHi(n-1)=a^(1/2)*(kk-X(n-1))-a^(3/2)/2*(kk-X(n-1))^2+1/2*a^(-1/2);
    PHi(1:n-2)=1/sqrt(a);

    price=W*beta0*PHi.'-S_0*beta1*PHi.';

end

end

