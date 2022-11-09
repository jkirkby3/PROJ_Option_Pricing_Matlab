function price = Price_GMDB_BSM_NoSurrender(T, M, gmdb_params, sigma, r, q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Analytical Price of Gauranteed Minimum death benefits with period fees and
%               NO early surrender under Black-Scholes model
%
% Author:      Justin Lars Kirkby
% References:  (1) Valuation and optimal surrender of variable annuities
%                   with guaranteed minimum benefits and periodic fees, 
%                   Kirkby and Aguilar 2022, Scandinavian Actuarial Journal
%
% Models Supported: Black-Scholes model
% Returns: price of contract
%
% ----------------------
% Contract/Model Params
% ----------------------
% T   = time remaining until maturity
% M   = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
%
% modelInputs
% ------------------------
% r     = interest rate
% q     = dividend yield
% sigma = constant Black-Scholes volatility
%
% gmdb_params
% ------------------------
% F_0       = initial fund value
% alpha_fee = period fee rate
% g         = floor on growth
% c         = cap on growth
% death_prob = probability of death (mortality) table, based on age. This
%           is NOT the conditional probability. To construct see function
%           make_mortality_table_pmf, set conditional=0
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% alph =  grid with is 2*alph
% N  = number of grid points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = T/M;
death_prob = gmdb_params.death_prob;
alpha = gmdb_params.alpha_fee;
drops = 1;  % drops due to the fee

F_0 = gmdb_params.F_0;

g = gmdb_params.g;
c = gmdb_params.c;  % Cap on growth 

price = 0;

for m = 1 : M
    t_m = m * dt;
    
    drops = drops * (1 - alpha);

    v_m = exp((r-q - sigma^2/2)*t_m) * drops; 
    a_m = g*t_m - log(v_m);
    b_m = c*t_m - log(v_m);
    
    sigtm = sigma^2*t_m;
    denom = sigma*sqrt(t_m);
    
    temp1 = normcdf((-a_m + sigtm)/denom) - normcdf((-b_m + sigtm)/denom);
    temp1 = v_m * exp(sigtm/2) * temp1;
    temp2 = exp(g*t_m)*normcdf(a_m/denom) + exp(c*t_m)*normcdf(-b_m/denom);
    
    price = price + death_prob(m) * exp(-r*t_m) * (temp1 + temp2);
end


price = price * F_0;
end

