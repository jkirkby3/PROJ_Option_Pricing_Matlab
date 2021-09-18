function modelInputs = getModelInput(model, dt, r, q, modelParams, T)
% model: Levy models (BSM, CGMY, NIG, MJD, Kou)
%        Affine models (Heston)
% r: interest rate
% q: dividend yield
% dt: time increment (this could be time to maturity, or time between monitoring dates)
% modelParams: dictionary of params for specific model
% NOTE: r,q can also be chosen to model forward or FX
% NOTE: the returned cumulants, c1, c2, c4, contain only contain dt in the
% case of Heston, but not the Levy models (this is so we can use this
% function in Exotic contexts)

if nargin < 6
   T = dt;  % Optional param allowing us to obtain an rnCHF(dt) and rnCHF(T), useful in some cases 
end

modelInputs = {};
modelInputs.dt = dt;
modelInputs.T = T;
modelInputs.r = r;
modelInputs.q = q;
modelInputs.model = model;

if model == 1 %BSM (Black Scholes Merton)
    %----------------------------------------------
    % Unpack the model parameters
    sigmaBSM = modelParams.sigmaBSM;
    %----------------------------------------------
    % Set the return object
    w = -.5*sigmaBSM^2;  % Convexity correction
    modelInputs.RNmu = r - q + w;   % DEF: Risk Neutral drift
    
    % Cumulants (useful for setting trunction range for density support / grids)
    modelInputs.c1 = modelInputs.RNmu;       % DEF: first cumulant
    modelInputs.c2 = sigmaBSM^2;             % DEF: second cumulant
    modelInputs.c4 = 0;                         % DEF: fourth cumulant
    
    % Charachteristic Functions / Levy Symbol
    modelInputs.rnCHF = @(u) cf_RN_BSM( u, r-q, dt, sigmaBSM );  % DEF: risk-neutral characteristic function at dt
    modelInputs.rnCHF_T = @(u) cf_RN_BSM( u, r-q, T, sigmaBSM );  % DEF: risk-neutral characteristic function at T
    modelInputs.rnSYMB = @(u) SYMB_RN_BSM(u, r-q, sigmaBSM  );  % DEF: risk neutral Levy symbol
    
elseif model == 2 %CGMY
    %----------------------------------------------
    % Unpack the model parameters
    C = modelParams.C; 
    G = modelParams.G; 
    MM = modelParams.MM; 
    Y = modelParams.Y; 
    %----------------------------------------------
    % Set the return object
    
    % Set the Risk-Neutral drift (based on interest/div rate and convexity correction)
    w = -C*gamma(-Y)*((MM-1)^Y - MM^Y + (G+1)^Y - G^Y);  % convexity correction
    modelInputs.RNmu = r - q + w;    % DEF: Risk Neutral drift
    
    % Cumulants (useful for setting trunction range for density support / grids)
    modelInputs.c1 = modelInputs.RNmu + C*gamma(1-Y)*(MM^(Y-1)-G^(Y-1));   % DEF: first cumulant
    modelInputs.c2 = C *gamma(2-Y)*(MM^(Y-2)+G^(Y-2));   % DEF: second cumulant
    modelInputs.c4 = C *gamma(4-Y)*(MM^(Y-4)+G^(Y-4));   % DEF: fourth cumulant
    
    % Charachteristic Functions / Levy Symbol
    modelInputs.rnCHF = @(u)cf_RN_CGMY(u,dt,r-q,C,G,MM,Y);  % DEF: risk-neutral characteristic function at dt
    modelInputs.rnCHF_T = @(u)cf_RN_CGMY(u,T,r-q,C,G,MM,Y);  % DEF: risk-neutral characteristic function at T
    modelInputs.rnSYMB = @(u) SYMB_RN_CGMY(u,r-q,C,G,MM,Y);  % DEF: risk neutral Levy symbol
    
elseif model == 3 %NIG  (Normal-Inverse-Gaussian)
    %----------------------------------------------
    alpha = modelParams.alpha;  
    beta = modelParams.beta;  
    delta = modelParams.delta;   
    %----------------------------------------------
    asq = alpha^2;  
    bsq = beta^2;  
    temp = sqrt(asq-bsq);
    
    % Set the Risk-Neutral drift (based on interest/div rate and convexity correction)
    w = delta*(sqrt(asq - (beta+1)^2)-temp);  % convexity correction
    modelInputs.RNmu = r - q + w;
    
    % Cumulants (useful for setting trunction range for density support / grids)
    modelInputs.c1 = modelInputs.RNmu + delta*beta/temp;
    modelInputs.c2 = delta*asq*(asq - bsq)^(-1.5);
    modelInputs.c4 = 3*delta*asq*(asq + 4*bsq)*(asq - bsq)^(-3.5);
    
    % Charachteristic Functions / Levy Symbol
    modelInputs.rnCHF = @(u) cf_RN_NIG( u,r-q,dt,alpha,beta,delta);
    modelInputs.rnCHF_T = @(u) cf_RN_NIG( u,r-q,T,alpha,beta,delta);
    modelInputs.rnSYMB = @(u) SYMB_RN_NIG(u,r-q,alpha,beta,delta);  % DEF: risk neutral Levy symbol
    
elseif model == 4 %MJD (Merton Jump Diffusion)
    %----------------------------------------------
    sigma = modelParams.sigma;  
    lam= modelParams.lam;  
    muj = modelParams.muj; 
    sigmaj = modelParams.sigmaj;
    %----------------------------------------------
    
    % Set the Risk-Neutral drift (based on interest/div rate and convexity correction)
    modelInputs.sig2 = .5*sigma^2; 
    w = -.5*sigma^2 - lam*(exp(muj + .5*sigmaj^2)-1);  % convexity correction
    modelInputs.RNmu = r - q + w;
    
    % Cumulants (useful for setting trunction range for density support / grids)
    modelInputs.c1 = modelInputs.RNmu + lam*muj;
    modelInputs.c2 = lam*(sigma^2/lam + muj^2 +sigmaj^2);
    modelInputs.c4 = lam*(muj^4 + 6*sigmaj^2*muj^2+3*sigmaj^4);
    
    % Charachteristic Functions / Levy Symbol
    modelInputs.rnCHF = @(u) cf_RN_MJD( u, r-q, dt, sigma, muj, sigmaj , lam);
    modelInputs.rnCHF_T = @(u) cf_RN_MJD( u, r-q, T, sigma, muj, sigmaj , lam);
    modelInputs.rnSYMB = @(u) SYMB_RN_MJD(u,r-q,sigma, muj, sigmaj , lam);  % DEF: risk neutral Levy symbol
    
elseif model == 5 %Kou's Double Expo
    %----------------------------------------------
    sigma = modelParams.sigma; 
    lam = modelParams.lam;
    p_up = modelParams.p_up; 
    eta1 = modelParams.eta1;
    eta2 = modelParams.eta2; 
    %----------------------------------------------
    
    % Set the Risk-Neutral drift (based on interest/div rate and convexity correction)
    w = -.5*sigma^2 - lam*(p_up*eta1/(eta1-1) + (1-p_up)*eta2/(eta2+1)-1);  % convexity correction
    modelInputs.RNmu = r - q + w;
    
    % Cumulants (useful for setting trunction range for density support / grids)
    modelInputs.c1 = modelInputs.RNmu + lam*p_up/eta1 + lam*(1-p_up)/eta2;
    modelInputs.c2 = sigma^2 + 2*lam*p_up/(eta1^2) + 2*lam*(1-p_up)/(eta2^2);
    modelInputs.c4 = 24*lam*(p_up/eta1^4 + (1-p_up)/eta2^4);
    
    % Charachteristic Functions / Levy Symbol
    modelInputs.rnCHF = @(u)cf_RN_KOU(u,dt,r-q,sigma,lam,p_up,eta1,eta2);
    modelInputs.rnCHF_T = @(u)cf_RN_KOU(u,T,r-q,sigma,lam,p_up,eta1,eta2);
    modelInputs.rnSYMB = @(u) SYMB_RN_Kou( u, r-q, sigma,lam,p_up,eta1,eta2);  % DEF: risk neutral Levy symbol
    
elseif model == 6  % Heston's model (Note: there are two supported parameter conventions here)
    %----------------------------------------------
    theta = modelParams.theta;
    rho = modelParams.rho; 
    
    if isfield(modelParams,'kappa')
        kappa = modelParams.kappa; 
    else
        kappa = modelParams.eta;
    end
    
    if isfield(modelParams,'v_0')
        v_0 = modelParams.v_0; 
    else
        v_0 = modelParams.v0;
    end
    
    if isfield(modelParams,'sigma_v')
        sigma_v = modelParams.sigma_v; 
    else
        sigma_v = modelParams.Sigmav;
    end

    %----------------------------------------------
    
    % Set the Risk-Neutral drift
    modelInputs.RNmu = r - q - .5*theta;
    
    % Cumulants (useful for setting trunction range for density support / grids)
    modelInputs.c1 = modelInputs.RNmu*dt + (1-exp(-kappa*dt))*(theta - v_0)/(2*kappa);
    modelInputs.c2 = 1/(8*kappa^3)*(sigma_v*dt*kappa*exp(-kappa*dt)*(v_0-theta)*(8*kappa*rho-4*sigma_v)...
        + kappa*rho*sigma_v*(1-exp(-kappa*dt))*(16*theta - 8*v_0)...
        +2*theta*kappa*dt*(-4*kappa*rho*sigma_v + sigma_v^2 + 4*kappa^2)...
        + sigma_v^2*((theta - 2*v_0)*exp(-2*kappa*dt) + theta*(6*exp(-kappa*dt)-7)+2*v_0)...
        + 8*kappa^2*(v_0-theta)*(1-exp(-kappa*dt)));
    modelInputs.c4 = 0;
    
    % Charachteristic Functions / Levy Symbol
    modelInputs.rnCHF = @(u)cf_RN_Heston(u,dt,r-q,v_0,theta,kappa,sigma_v,rho);
    modelInputs.rnCHF_T = @(u)cf_RN_Heston(u,T,r-q,v_0,theta,kappa,sigma_v,rho);
    % NOTE: no rnSYMB for this model, as we have no current use for it
    
elseif model == 7 %KoBoL Model  
    %----------------------------------------------
    % Unpack the model parameters
    c = modelParams.c; 
    lam_p = modelParams.lam_p; 
    lam_m = modelParams.lam_m; 
    nu = modelParams.nu; 
    %----------------------------------------------
    
    % NOTE: params have been
    % written in correspondence with CGMY, which is a subclass of KoBoL
    C = c; MM = lam_p; G = -lam_m; Y = nu;
    
    % Set the Risk-Neutral drift (based on interest/div rate and convexity correction)
    w = -C*gamma(-Y)*((MM-1)^Y - MM^Y + (G+1)^Y - G^Y);  % convexity correction
    modelInputs.RNmu = r - q + w;    % DEF: Risk Neutral drift
    
    % Cumulants (useful for setting trunction range for density support / grids)
    modelInputs.c1 = modelInputs.RNmu + C*gamma(1-Y)*(MM^(Y-1)-G^(Y-1));   % DEF: first cumulant
    modelInputs.c2 = C *gamma(2-Y)*(MM^(Y-2)+G^(Y-2));   % DEF: second cumulant
    modelInputs.c4 = C *gamma(4-Y)*(MM^(Y-4)+G^(Y-4));   % DEF: fourth cumulant
    
    % Charachteristic Functions / Levy Symbol
    modelInputs.rnCHF = @(u)cf_RN_KoBoL(u,dt,r-q,c,lam_p,lam_m,nu);  % DEF: risk-neutral characteristic function at dt
    modelInputs.rnCHF_T = @(u)cf_RN_KoBoL(u,T,r-q,c,lam_p,lam_m,nu);  % DEF: risk-neutral characteristic function at T
    modelInputs.rnSYMB = @(u) SYMB_RN_KoBoL(u,r-q,c,lam_p,lam_m,nu);  % DEF: risk neutral Levy symbol
    
elseif model == 8 %Variance Gamma
    %----------------------------------------------
    % Unpack the model parameters
    sigma = modelParams.sigma;  
    theta = modelParams.theta;  
    nu = modelParams.nu; 
    %----------------------------------------------
    
    % Set the Risk-Neutral drift (based on interest/div rate and convexity correction)
    sig2 = .5*sigma^2; 
    w = log(1 - theta*nu - sig2*nu)/nu; % convexity correction
    modelInputs.RNmu = r - q + w;
    
    % Cumulants (useful for setting trunction range for density support / grids)
    modelInputs.c1 = modelInputs.RNmu + theta;
    modelInputs.c2 = (sigma*sigma + nu*theta*theta);
    modelInputs.c4 = 3*(sigma^4*nu + 2*theta^4*nu^3 + 4*sigma^2*theta^2*nu^2);
    
    % Charachteristic Functions / Levy Symbol
    modelInputs.rnCHF = @(u) cf_RN_VG( u, r-q, dt, sigma, nu, theta);
    modelInputs.rnCHF_T = @(u) cf_RN_VG( u, r-q, T, sigma, nu, theta);
    modelInputs.rnSYMB = @(u) SYMB_RN_VG(u, r-q, sigma, nu, theta);  % DEF: risk neutral Levy symbol
    
elseif model == 9 %Bilateral Gamma
    %----------------------------------------------
    % Unpack the model parameters
    alpha_p = modelParams.alpha_p;
    lam_p = modelParams.lam_p;
    alpha_m = modelParams.alpha_m;
    lam_m = modelParams.lam_m;
    %----------------------------------------------
    
    % Set the Risk-Neutral drift (based on interest/div rate and convexity correction)
    m1 = alpha_p / lam_p - alpha_m / lam_m;
    w = -log((lam_p/(lam_p -1))^alpha_p*(lam_m/(lam_m +1))^alpha_m); % convexity correction

    modelInputs.RNmu = r - q + w;
    
    cumulants = @(n) factorial(n-1)*(alpha_p/lam_p^n + (-1)^n*alpha_m/lam_m^n);
    
    % Cumulants (useful for setting trunction range for density support / grids)
    modelInputs.c1 = modelInputs.RNmu + m1;
    modelInputs.c2 = cumulants(2);
    modelInputs.c4 = cumulants(4);
    
    % Charachteristic Functions / Levy Symbol
    modelInputs.rnCHF = @(u) cf_RN_BilateralGamma( u, r-q, dt, alpha_p, lam_p, alpha_m, lam_m);
    modelInputs.rnCHF_T = @(u) cf_RN_BilateralGamma( u, r-q, T, alpha_p, lam_p, alpha_m, lam_m);
    modelInputs.rnSYMB = @(u) SYMB_RN_BilateralGamma(u, r-q, alpha_p, lam_p, alpha_m, lam_m);  % DEF: risk neutral Levy symbol
end



end


