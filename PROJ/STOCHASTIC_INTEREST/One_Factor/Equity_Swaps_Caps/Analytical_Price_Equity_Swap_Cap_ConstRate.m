function [ price ] = Analytical_Price_Equity_Swap_Cap_ConstRate(r, M, dt, contract_param, sigma )
% Price the equity swap and cap/floor under constant interest rate assumptino

contract_type = contract_param.contract_type; % 1 = swap, 2 = cap/floor
K = contract_param.K;
c = contract_param.c;
T = dt*M;

gamma = (1 + r*dt + c);
mult = (1 - exp(-r*T))/(exp(r*dt) - 1);


if contract_type == 1  % Swap
    price = K*(exp(r*dt) - gamma)*mult;
    
elseif contract_type == 2  % cap/floor
    
    C = contract_param.C; % cap
    F = contract_param.F; % floor
    
    mu = (r - sigma*sigma/2)*dt;
    sigdt = sigma*sqrt(dt);
    a = (log(max(0, F + gamma)) - mu)/sigdt;
    b = (log(C + gamma) - mu)/sigdt;
    t1 = normcdf(b - sigdt) - normcdf(a - sigdt);
    t2 = normcdf(b) - normcdf(a);
    H = F*normcdf(a) + exp(r*dt)*t1 - gamma*t2 + C*(1-normcdf(b));
    price = K*H*mult;
else
    price = -123456789;
end

end

