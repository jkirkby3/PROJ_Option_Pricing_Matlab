function price = BinomialLattice_BlackScholes_func(S_0, K, r, T, sigma, M, call, american)
% European and American/Bermudan Option pricer 
% Model: Black Scholes Merton
% Author: Justin Kirkby

dt = T/M;
u = exp(sigma*sqrt(dt));
d = 1/u;
p = (exp(r*dt)-d)/(u-d);
discount = exp(-r*dt);
p_u = discount*p;
p_d = discount*(1-p);

% Determine Payoff Function
if call == 1
    payoff = @(C, S)  max( C, S - K);  % C is the continuation value
else
    payoff = @(C, S)  max( C, K - S);  % C is the continuation value
end

Svals = zeros(2*M +1, 1);
Svals(M+1) = S_0;
for i = 1:M
    Svals(M+1+i) = u*Svals(M+i);
    Svals(M+1-i) = d*Svals(M+2-i);
end

% Initialize Terminal Payoff
Pvals = zeros(2*M + 1, 1);
for i=1:2:2*M +1
    Pvals(i) = payoff(0, Svals(i));
end

% Calculate Price Recursively
if american == 1
    for tau = 1:M
        for i = (tau+1):2:(2*M +1-tau)
            cont = p_u*Pvals(i+1) + p_d*Pvals(i-1);
            Pvals(i) = payoff(cont, Svals(i));
        end
    end
else % European
    for tau = 1:M
        for i = (tau+1):2:(2*M +1-tau)
            Pvals(i) = p_u*Pvals(i+1) + p_d*Pvals(i-1);
        end
    end
end
    
price = Pvals(M+1);

end

