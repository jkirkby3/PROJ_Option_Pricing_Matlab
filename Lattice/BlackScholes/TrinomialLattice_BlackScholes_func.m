function price = TrinomialLattice_BlackScholes_func(S_0, K, r, T, sigma, M, call, american)
% European and American/Bermudan Option pricer 
% Model: Black Scholes Merton
% Author: Justin Kirkby

dt = T/M;
dx = sigma*sqrt(dt);  % this could be chosen differently
s2 = r - 0.5*sigma^2;
discount = exp(-r*dt);
cons1 = (sigma^2*dt + (s2*dt)^2) / dx^2;
cons2 = s2*dt/dx;

p_u = discount * 0.5 * (cons1 + cons2);
p_m = discount * (1 - cons1);
p_d = discount * 0.5 * (cons1 - cons2);


% Determine Payoff Function
if call == 1
    payoff = @(C, S)  max( C, S - K);  % C is the continuation value
else
    payoff = @(C, S)  max( C, K - S);  % C is the continuation value
end

edx = exp(dx);

% Initialize terminal S grid values
Svals = zeros(2*M +1, 1);
Svals(1) = S_0 * exp(-M*dx);
for i = 2:2*M + 1
    Svals(i) = edx*Svals(i-1);
end

% Initialize Terminal Payoff
Pvals = zeros(2*M + 1, 2);
tau = mod(M, 2) + 1;
for i=1:2*M +1
    Pvals(i, tau) = payoff(0, Svals(i));
end

% Calculate Price Recursively
if american == 1
    for tau = M-1:-1:0
        k = mod(tau,2) + 1;
        k1 = mod(tau + 1, 2) + 1;
        for i = (M - tau + 1): (M + tau + 1)
            cont = p_d * Pvals(i-1,k1) + p_m * Pvals(i, k1) + p_u * Pvals(i+1, k1);
            Pvals(i, k) = payoff(cont, Svals(i));
        end
    end
else % European
    for tau = M-1:-1:0
        k = mod(tau,2) + 1;
        k1 = mod(tau + 1, 2) + 1;
        for i = (M - tau + 1): (M + tau + 1)
            Pvals(i, k) = p_d * Pvals(i-1,k1) + p_m * Pvals(i, k1) + p_u * Pvals(i+1, k1);
        end
    end
end
    
price = Pvals(M+1,1);

end

