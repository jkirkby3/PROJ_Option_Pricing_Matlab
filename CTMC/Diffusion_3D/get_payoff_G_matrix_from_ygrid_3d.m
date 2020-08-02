function G = get_payoff_G_matrix_from_ygrid_3d( y_1, y_2, y_3, S_0s, sigmas, R, contractParams)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

payoff_type = contractParams.payoff_type;

rho12 = R(1,2);
rho23 = R(2,3);
rho13 = R(1,3);

gamma = (rho12*rho13 - rho23)/(1 - rho12^2);

if payoff_type == 1 || payoff_type == 2   % G = S_1, or G=S_2, or G=S_3
    dim = contractParams.dim;
    if dim == 1
        payoff = @(y1,y2,y3)S_0s(1)*exp(sigmas(1)*y1);
    elseif dim == 2
        payoff = @(y1,y2,y3)S_0s(2)*exp(sigmas(2)*(y2 + rho12*y1));
    else
        payoff = @(y1,y2,y3)S_0s(3)*exp(sigmas(3)*(rho13*y1 -gamma*y2 + y3));
    end
    

elseif payoff_type == 5  % Geometric Basket Call / Put:  G = (sqrt(S_1) * sqrt(S_2) - K)^+ and
    K = contractParams.K;
    if contractParams.call == 1
        payoff = @(y1,y2,y3) max(0, (S_0s(1)*exp(sigmas(1)*y1) * S_0s(2)*exp(sigmas(2)*(y2 + rho12*y1))* S_0s(3)*exp(sigmas(3)*(rho13*y1 -gamma*y2 + y3)))^(1/3) - K);
    else
        payoff = @(y1,y2,y3) max(0, K - (S_0s(1)*exp(sigmas(1)*y1) * S_0s(2)*exp(sigmas(2)*(y2 + rho12*y1))* S_0s(3)*exp(sigmas(3)*(rho13*y1 -gamma*y2 + y3)))^(1/3)); 
    end
    
elseif payoff_type == 6  % Arithmetic Basket Call / Put:  G = (sqrt(S_1) * sqrt(S_2) - K)^+
    K = contractParams.K;
    if contractParams.call == 1
        payoff = @(y1,y2,y3) max(0, (1/3)*(S_0s(1)*exp(sigmas(1)*y1) + S_0s(2)*exp(sigmas(2)*(y2 + rho12*y1)) + S_0s(3)*exp(sigmas(3)*(rho13*y1 -gamma*y2 + y3))) - K);
    else
        payoff = @(y1,y2,y3) max(0, K - (1/3)*(S_0s(1)*exp(sigmas(1)*y1) + S_0s(2)*exp(sigmas(2)*(y2 + rho12*y1)) + S_0s(3)*exp(sigmas(3)*(rho13*y1 -gamma*y2 + y3)))); 
    end
    

    
end
    
m_0 = length(y_1);
G = zeros(m_0, m_0, m_0);

for i=1:m_0
    for j=1:m_0
        for k=1:m_0
            G(i,j,k) = payoff(y_1(i), y_2(j), y_3(k));
        end
    end
end

end

