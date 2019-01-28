function EXP_A = get_SV_matrix_exponential( Q, dt, xi, v1, v2, fv, psi_J, m_0, N )

Qt = dt*Q';
EXP_A = ones(m_0,m_0,N);
for j = 1:N  
    EXP_A(:,:,j) = Qt + diag(v1*xi(j) - v2*xi(j)^2 + dt*psi_J(xi(j)));  %This incorporates the dt
end

%%% Define Toeplitz matrix (note the special structure:
%%% .. Before exponentiation it is: Toeplitz + skew-symmetric + zero diagonal
Lambda_dxi = zeros(m_0,m_0);

for k=1:m_0
    for j = k+1:m_0
        Lambda_dxi(k,j) = fv(k) - fv(j);
        Lambda_dxi(j,k) = -Lambda_dxi(k,j);
    end
end

Lambda_dxi = exp(Lambda_dxi); %note: not matrix exponential

EXP_A(:,:,1) = expm(EXP_A(:,:,1));  %matrix exponential 
%%% (note: for n=1, Lambda is matrix of ones, so we start next loop at n=2)

Lambda = Lambda_dxi; %initialize Lambda for n = 2

for j = 2:N 
    EXP_A(:,:,j) = expm(EXP_A(:,:,j)).*Lambda;  %matrix exponential 
    Lambda = Lambda.*Lambda_dxi; %update Lambda
end

end

