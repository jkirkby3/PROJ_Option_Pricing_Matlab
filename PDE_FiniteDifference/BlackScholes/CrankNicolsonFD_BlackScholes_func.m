function price = CrankNicolsonFD_BlackScholes_func( S_0, K, r, T, sigma, call, dS, dt, Smax, Smin)
% Description: Crank-Nicolson PDE Finite Difference method to price European Option in Black Scholes Model
% Author: Justin Kirkby
M = round((Smax - Smin)/dS); % grid points
N = round(T/dt); % time steps

dS = (Smax - Smin) / M; % readjust
dt = T / N;  % readjust

vals = zeros(M+1, N+1);
vS = linspace(Smin, Smax, M+1)';
vI = vS / dS;
vJ = 0:N;

% Boundary Conditions
if call == 1  % call option
    vals(:, N+1) = max(vS - K, 0);  
    vals(1, :) = 0;
    vals(M+1, :) = Smax - K*exp(-r*dt*(N - vJ)); 
else % put option
    vals(:, N+1) = max(K - vS, 0);  
    vals(1, :) = K*exp(-r*dt*(N - vJ));
    vals(M+1, :) = 0; 
end

% Tridiagonal Coefficients
a = 0.25 * dt * (sigma^2 *(vI.^2) - r*vI);
b = -dt*0.5*(sigma^2*(vI.^2) + r);
c = 0.25*dt*(sigma^2*(vI.^2) + r*vI);

M1 = -diag(a(3:M), -1) + diag(1 - b(2:M)) - diag(c(2:M-1),1);
[L,U] = lu(M1);
M2 = diag(a(3:M), -1) + diag(1 + b(2:M)) + diag(c(2:M-1),1);

% Solve systems (backward in time)
for j = N:-1:1
    vals(2:M,j) = U \ (L \ (M2*vals(2:M,j+1)));
end

price = interp1(vS, vals(:,1), S_0);

end

