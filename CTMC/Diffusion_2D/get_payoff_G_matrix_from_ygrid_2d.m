function G = get_payoff_G_matrix_from_ygrid_2d( y_1, y_2, S_0s, sigmas, rho, contractParams)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

payoff_type = contractParams.payoff_type;

if payoff_type == 1  % G = S_1
    payoff = @(y1,y2)S_0s(1)*exp(sigmas(1)*y1);  
    
elseif payoff_type == 2  % G = S_2
    payoff = @(y1,y2)S_0s(2)*exp(sigmas(2)*(y2 + rho*y1));
    
elseif payoff_type == 3  % Exchange:  G = (S_1 - S_2)^+
    payoff = @(y1,y2) max(0, S_0s(1)*exp(sigmas(1)*y1) - S_0s(2)*exp(sigmas(2)*(y2 + rho*y1)));
    
elseif payoff_type == 4  % Spread:  G = (S_1 - S_2 - K)^+
    K = contractParams.K;
    payoff = @(y1,y2) max(0, S_0s(1)*exp(sigmas(1)*y1) - S_0s(2)*exp(sigmas(2)*(y2 + rho*y1)) - K);
    
elseif payoff_type == 5  % Geometric Basket Call / Put:  G = (sqrt(S_1) * sqrt(S_2) - K)^+ and
    K = contractParams.K;
    if contractParams.call == 1
        payoff = @(y1,y2) max(0, sqrt(S_0s(1)*exp(sigmas(1)*y1)) * sqrt(S_0s(2)*exp(sigmas(2)*(y2 + rho*y1))) - K);
    else
        payoff = @(y1,y2) max(0, K - sqrt(S_0s(1)*exp(sigmas(1)*y1)) * sqrt(S_0s(2)*exp(sigmas(2)*(y2 + rho*y1)))); 
    end
    
elseif payoff_type == 6  % Arithmetic Basket Call / Put:  G = (sqrt(S_1) * sqrt(S_2) - K)^+
    K = contractParams.K;
    if contractParams.call == 1
        payoff = @(y1,y2) max(0, 0.5*S_0s(1)*exp(sigmas(1)*y1) + 0.5*S_0s(2)*exp(sigmas(2)*(y2 + rho*y1)) - K);
    else
        payoff = @(y1,y2) max(0, K - 0.5*S_0s(1)*exp(sigmas(1)*y1) - 0.5*S_0s(2)*exp(sigmas(2)*(y2 + rho*y1))); 
    end
    
elseif payoff_type == 7  % Call-on-Max and Put-on-Min
    K = contractParams.K;
    if contractParams.call == 1
        payoff = @(y1,y2) max(0, max(S_0s(1)*exp(sigmas(1)*y1), S_0s(2)*exp(sigmas(2)*(y2 + rho*y1))) - K);
    else
        payoff = @(y1,y2) max(0, K - min(S_0s(1)*exp(sigmas(1)*y1), S_0s(2)*exp(sigmas(2)*(y2 + rho*y1)))); 
    end
elseif payoff_type == 8  % Call/put on Just S_2
    K = contractParams.K;
    if contractParams.call == 1
        payoff = @(y1,y2) max(0, S_0s(2)*exp(sigmas(2)*(y2 + rho*y1)) - K);
    else
        payoff = @(y1,y2) max(0, K - S_0s(2)*exp(sigmas(2)*(y2 + rho*y1))); 
    end    
elseif payoff_type == 9  % Best-of / worst of
    if contractParams.best == 1
        payoff = @(y1,y2) max(S_0s(1)*exp(sigmas(1)*y1), S_0s(2)*exp(sigmas(2)*(y2 + rho*y1)));
    else  % worst of
        payoff = @(y1,y2) min(S_0s(1)*exp(sigmas(1)*y1), S_0s(2)*exp(sigmas(2)*(y2 + rho*y1)));
    end
end
    
m_0 = length(y_1);
G = zeros(m_0, m_0);

for i=1:m_0
    for j=1:m_0
        G(i,j) = payoff(y_1(i), y_2(j));
    end
end

end

