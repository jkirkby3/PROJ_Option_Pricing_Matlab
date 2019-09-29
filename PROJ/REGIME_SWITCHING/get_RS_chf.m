function chf = get_RS_chf( Q, dt, xi, drifts, vols, initial_state)
% Risk-Neutral Chf of Regime Switching model with given initial_state
N = length(xi);

Qt = dt*Q;
chf = zeros(size(xi));

vols = 0.5*vols.^2;
drifts = dt*(drifts - vols)*1i;

% TODO: vectorize
for j = 1:N  
    temp = expm(Qt + diag(drifts*xi(j) - vols*xi(j)^2));
    chf(j) = sum(temp(:, initial_state));
end

end