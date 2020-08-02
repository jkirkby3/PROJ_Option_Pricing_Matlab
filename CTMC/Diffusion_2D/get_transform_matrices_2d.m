function [ L, D, C, Cinv ] = get_transform_matrices_2d( R, method )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if method == 1
    % LDL Without the permutation
    [L, D] = ldl(R);
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
    rho = R(1,2);
    L = [1 0; rho 1];
    D = diag([1, 1 - rho^2]);
    C = [1 0; -rho 1]; 
    % C = inv(L);
    Cinv = L;
end


end

