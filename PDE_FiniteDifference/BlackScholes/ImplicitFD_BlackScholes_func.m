function price = ImplicitFD_BlackScholes_func( S_0, K, r, T, sigma, call, dS, dt, Smax, Smin)
% Description: Fully Implicit PDE Finite Difference method to price European Option in Black Scholes Model
% Author: Justin Kirkby
M = round((Smax - Smin)/dS); % grid points
N = round(T/dt); % time steps

% dS = (Smax - Smin) / M; % readjust
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
a = 0.5 * (r * dt * vI - sigma^2*dt*(vI.^2));
b = 1 + sigma^2*dt*(vI.^2) + r*dt;
c = - 0.5 * (r * dt * vI + sigma^2*dt*(vI.^2));

cvec = diag(a(3:M), -1) + diag(b(2:M)) + diag(c(2:M-1),1);
[L,U] = lu(cvec);

% Solve systems (backward in time)
z = zeros(M-1,1);
for j = N:-1:1
    z(1) = -a(2) * vals(1,j);
    vals(2:M,j) = U \ (L \ (vals(2:M,j+1) + z));
end

price = interp1(vS, vals(:,1), S_0);

end

