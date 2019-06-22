function price = ExplicitFD_BlackScholes_func( S_0, K, r, T, sigma, call, dS, dt, Smax, Smin)
% Description: Explicit PDE Finite Difference method to price European Option in Black Scholes Model
% Author: Justin Kirkby
% NOTE: to ensure stability, we can choose dS, and then dt = dS^2 / Smax^2 / sigma^2;
M = round((Smax - Smin)/dS); % grid points
N = round(T/dt); % time steps

% dS = Smax / M; % readjust
dt = T / N;  % readjust

vals = zeros(M+1, N+1);
vS = linspace(Smin, Smax, M+1)';
vI = 0:M;
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
a = 0.5 * dt * (sigma^2*vI - r).*vI;
b = 1 - dt * (sigma^2*vI.^2 + r);
c = 0.5 * dt * (sigma^2*vI + r).*vI;

% Solve (backward)
for j = N:-1:1
   for i=2:M
      vals(i,j) = a(i)*vals(i-1,j+1) + b(i)*vals(i,j+1) + c(i)*vals(i+1,j+1); 
   end
end

price = interp1(vS, vals(:,1), S_0);

end

