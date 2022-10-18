function [ price, r_mean, r_var, mart] = Price_Bond_CTMC_OneFactor( numeric_param, T, model, modparam, do_moments)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Price Bonds under 1 Factor short rate model
% Models Supported: 1 Factor short rate model
% Returns: price of bond
% Author: Justin Lars Kirkby
%
% References:  (1) 
%
% ----------------------
% Contract Params 
% ----------------------
% T  : time remaining to maturity (T = 2 is two years, T = .5 is half a year)
%
% ----------------------
% Model Params 
% ----------------------
% model:
%        1 = HESTON:      Sigmar, r0, rho, eta, theta
%        2 = VASICEK:     Sigmar, r0, rho, eta, theta
%
% modparam: contains all necessary params for the specific model (see below during assingment which ones are needed)
%
% ----------------------
% Numerical Params 
% ----------------------
% numeric_parm: container of numerical params
%   m_0: number of states to approximate the Heston model with
%   gamma: var grid width parameter, grid is +/- gamma*stddev(variance process)
%   gridMethod: which type of var grid to use (typcially use 4)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_0 = numeric_param.m_0;
gridMethod = numeric_param.gridMethod;
gamma = numeric_param.gamma;
gridMultParam = numeric_param.gridMultParam;

if nargin < 5
    do_moments = 1;
end

t = T/2;
[lx, r0, ux] = get_rate_grid_boundaries( model, modparam, t, gamma);

[ mu_func,  sig_func,  zeta_func, varrho_func] = get_short_rate_grid_diffusion_funcs( model,  modparam);
boundaryMethod = 1;
center = r0; %this is where grid clusters... we can experiment with other choices.. 

[Q, rs, c_index] = Q_Matrix(m_0, mu_func,sig_func,lx,ux,gridMethod,center, gridMultParam, boundaryMethod);

prices = expm(T*(Q - diag(rs)))*ones(m_0,1);
price = prices(c_index);
% MM = expm(T*(Q - diag(rs))); 
% price = sum(MM(c_index,:));

if do_moments == 1
    P = expm(T*Q);
    means = P*rs.';
    r_mean = means(c_index);

    % plot(v, means)

    vars = P*(rs.*rs).' - means.*means;
    r_var = vars(c_index);

    % plot(v, vars)

    dt = T;
    Gamm = expm(dt*(Q - diag(varrho_func(rs))));
    rho = modparam.rho;
    Sigma_s = modparam.Sigma_s;
    kk = c_index;
    mart = 0;

    c1 = (rho*Sigma_s)^2/2*dt;
    for l=1:m_0
       mart = mart + exp(zeta_func(rs(kk), rs(l)) - c1)*Gamm(kk,l);
    end
else
    r_mean = 0; r_var = 0; mart = 0;
end

end

