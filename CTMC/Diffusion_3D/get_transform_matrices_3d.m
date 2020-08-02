function [ L, D, C, Cinv ] = get_transform_matrices_3d( R, method )

if method == 1
    % LDL Without the permutation
    [L,D] = ldl(R);
    C = inv(L); 
    Cinv = L;

elseif method == 2
    % Eigen Decomp
    [V,D] = eig(R);
    C = V.';
    Cinv = V;
    
elseif method == 3
    % Cholesky Decomp (Special case of LDL, D=I)
    L = chol(R);  % NOTE: this is the UPPER triangular transform
    L = L.';
    C = inv(L);
    Cinv = L;
    D = eye(2,2);

elseif method == 4  %%% ANALYTICAL
   rho12 = R(1,2);
   rho23 = R(2,3);
   rho13 = R(1,3);
    
   gamma = (rho12*rho13 - rho23)/(1 - rho12^2);
   L = [1 0 0; rho12 1 0; rho13 -gamma 1];
   C = [1 0 0; -rho12 1 0; (-gamma*rho12 - rho13) gamma 1];
   D = diag([1, 1 - rho12^2, det(R)/(1-rho12^2)]);
   Cinv = L;  

end

end

