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
if model == 1 %BSM (Black Scholes Merton)
    %----------------------------------------------
    % Unpack the model parameters
    sigmaBSM = modelParams.sigmaBSM;
    %----------------------------------------------
    % Set the return object
    modelInputs.RNmu = r - q - .5*sigmaBSM^2;   % DEF: Risk Neutral drift
    modelInputs.c1 = modelInputs.RNmu;       % DEF: first cumulant
    modelInputs.c2 = sigmaBSM^2;             % DEF: second cumulant
    modelInputs.c4 = 0;                         % DEF: fourth cumulant
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
    modelInputs.RNmu = r - q - C*gamma(-Y)*((MM-1)^Y - MM^Y + (G+1)^Y - G^Y);    % DEF: Risk Neutral drift
    modelInputs.c1 = modelInputs.RNmu + C*gamma(1-Y)*(MM^(Y-1)-G^(Y-1));   % DEF: first cumulant
    modelInputs.c2 = C *gamma(2-Y)*(MM^(Y-2)+G^(Y-2));   % DEF: second cumulant
    modelInputs.c4 = C *gamma(4-Y)*(MM^(Y-4)+G^(Y-4));   % DEF: fourth cumulant
    modelInputs.rnCHF = @(u)cf_RN_CGMY(u,dt,r-q,C,G,MM,Y);  % DEF: risk-neutral characteristic function at dt
    modelInputs.rnCHF_T = @(u)cf_RN_CGMY(u,T,r-q,C,G,MM,Y);  % DEF: risk-neutral characteristic function at T
    modelInputs.rnSYMB = @(u) SYMB_RN_CGMY(u,r-q,C,G,MM,Y);  % DEF: risk neutral Levy symbol
    
elseif model == 3 %NIG
    %----------------------------------------------
    alpha = modelParams.alpha;  
    beta = modelParams.beta;  
    delta = modelParams.delta;   
    %----------------------------------------------
    asq = alpha^2;  
    bsq = beta^2;  
    temp = sqrt(asq-bsq);
    
    modelInputs.RNmu = r - q + delta*(sqrt(asq - (beta+1)^2)-temp);
    modelInputs.c1 = modelInputs.RNmu + delta*beta/temp;
    modelInputs.c2 = delta*asq*(asq - bsq)^(-1.5);
    modelInputs.c4 = 3*delta*asq*(asq + 4*bsq)*(asq - bsq)^(-3.5);
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
    modelInputs.sig2 = .5*sigma^2; 
    modelInputs.RNmu = r - q - .5*sigma^2-lam*(exp(muj - .5*sigmaj^2)-1);
    modelInputs.c1 = modelInputs.RNmu + lam*muj;
    modelInputs.c2 = lam*(sigma^2/lam + muj^2 +sigmaj^2);
    modelInputs.c4 = lam*(muj^4 + 6*sigmaj^2*muj^2+3*sigmaj^4);
    modelInputs.rnCHF = @(u) cf_RN_MJD( u, r-q, dt, sigma, muj, sigmaj , lam);
    modelInputs.rnCHF_T = @(u) cf_RN_MJD( u, r-q, T, sigma, muj, sigmaj , lam);
    modelInputs.rnSYMB = @(u) SYMB_RN_MJD(u,r-q,sigma, muj, sigmaj , lam);  % DEF: risk neutral Levy symbol
    
elseif model == 5 %Kou Double Expo
    %----------------------------------------------
    sigma = modelParams.sigma; 
    lam = modelParams.lam;
    p_up = modelParams.p_up; 
    eta1 = modelParams.eta1;
    eta2 = modelParams.eta2; 
    %----------------------------------------------
    
    modelInputs.RNmu = r-q - .5*sigma^2 - lam*(p_up*eta1/(eta1-1) + (1-p_up)*eta2/(eta2+1)-1);
    modelInputs.c1 = modelInputs.RNmu + lam*p_up/eta1 + lam*(1-p_up)/eta2;
    modelInputs.c2 = sigma^2 + 2*lam*p_up/(eta1^2) + 2*lam*(1-p_up)/(eta2^2);
    modelInputs.c4 = 24*lam*(p_up/eta1^4 + (1-p_up)/eta2^4);
    modelInputs.rnCHF = @(u)cf_RN_KOU(u,dt,r-q,sigma,lam,p_up,eta1,eta2);
    modelInputs.rnCHF_T = @(u)cf_RN_KOU(u,T,r-q,sigma,lam,p_up,eta1,eta2);
    modelInputs.rnSYMB = @(u) SYMB_RN_Kou( u, r-q, sigma,lam,p_up,eta1,eta2);  % DEF: risk neutral Levy symbol
    
elseif model == 6  % Heston's model
    %----------------------------------------------
    v_0 = modelParams.v_0; 
    theta = modelParams.theta;
    kappa = modelParams.kappa; 
    sigma_v = modelParams.sigma_v;
    rho = modelParams.rho; 
    %----------------------------------------------
    modelInputs.RNmu = r-q -.5*theta;
    modelInputs.c1 = modelInputs.RNmu*dt + (1-exp(-kappa*dt))*(theta - v_0)/(2*kappa);
    modelInputs.c2 = 1/(8*kappa^3)*(sigma_v*dt*kappa*exp(-kappa*dt)*(v_0-theta)*(8*kappa*rho-4*sigma_v)...
        + kappa*rho*sigma_v*(1-exp(-kappa*dt))*(16*theta - 8*v_0)...
        +2*theta*kappa*dt*(-4*kappa*rho*sigma_v + sigma_v^2 + 4*kappa^2)...
        + sigma_v^2*((theta - 2*v_0)*exp(-2*kappa*dt) + theta*(6*exp(-kappa*dt)-7)+2*v_0)...
        + 8*kappa^2*(v_0-theta)*(1-exp(-kappa*dt)));
    modelInputs.c4 = 0;
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
    
    % Set the return object
    modelInputs.RNmu = r - q - C*gamma(-Y)*((MM-1)^Y - MM^Y + (G+1)^Y - G^Y);    % DEF: Risk Neutral drift
    modelInputs.c1 = modelInputs.RNmu + C*gamma(1-Y)*(MM^(Y-1)-G^(Y-1));   % DEF: first cumulant
    modelInputs.c2 = C *gamma(2-Y)*(MM^(Y-2)+G^(Y-2));   % DEF: second cumulant
    modelInputs.c4 = C *gamma(4-Y)*(MM^(Y-4)+G^(Y-4));   % DEF: fourth cumulant
    modelInputs.rnCHF = @(u)cf_RN_KoBoL(u,dt,r-q,c,lam_p,lam_m,nu);  % DEF: risk-neutral characteristic function at dt
    modelInputs.rnCHF_T = @(u)cf_RN_KoBoL(u,T,r-q,c,lam_p,lam_m,nu);  % DEF: risk-neutral characteristic function at T
    modelInputs.rnSYMB = @(u) SYMB_RN_KoBoL(u,r-q,c,lam_p,lam_m,nu);  % DEF: risk neutral Levy symbol
end



end


