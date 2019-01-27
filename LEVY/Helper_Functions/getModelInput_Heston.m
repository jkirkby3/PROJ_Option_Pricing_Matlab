function modelInputs = getModelInput_Heston(dt, r, q, modelParams)
% model: Levy models (BSM, CGMY, NIG, MJD, Kou)
%        Affine models (Heston)
% r: interest rate
% q: dividend yield
% dt: time increment (this could be time to maturity, or time between monitoring dates)
% modelParams: dictionary of params for specific model
% NOTE: r,q can also be chosen to model forward or FX

modelInputs = {};

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
% NOTE: no rnSYMB for this model, as we have no current use for it


end


