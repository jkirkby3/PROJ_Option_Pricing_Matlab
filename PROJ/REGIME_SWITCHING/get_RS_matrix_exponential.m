function EXP_A = get_RS_matrix_exponential( Q, dt, xi, drifts, vols, psi_J)
N = length(xi);
m_0 = length(drifts);

Qt = dt*Q';
EXP_A = ones(m_0,m_0,N);

drifts = dt*drifts;
vols = (0.5*dt)*vols.^2;

if nargin < 8
    for j = 1:N  
        EXP_A(:,:,j) = expm(Qt + diag(drifts*xi(j) - vols*xi(j)^2));  %This incorporates the dt
    end
else
    for j = 1:N  
        EXP_A(:,:,j) = expm(Qt + diag(drifts*xi(j) - vols*xi(j)^2 + dt*psi_J(xi(j))));  %This incorporates the dt
    end
end


end