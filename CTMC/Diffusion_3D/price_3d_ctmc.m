function [vals, c_index_1, c_index_2, c_index_3, y_1, y_2, y_3] = price_3d_ctmc( S_0s, T, r, R, sigmas, qs, params, contractParams, M)

if nargin < 9
    M = 1; % M is only needed for 
end
dt = T/M;

contract = contractParams.contract;

if contract == 1  % European
    dt = 1; M = 1;
end

method = 4;
num_devs = params.num_devs;
m_0 = params.m_0;
GridMultParam = params.GridMultParam;
gridMethod = params.gridMethod;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drifts = r - qs;

[ L, D, C, Cinv ] = get_transform_matrices_3d( R, method );

% Now Define New Uncorrelated System  (the dc underscore)
[drift_dc, sigma_dc ] = decorrelate(sigmas, drifts, C, D );

[Ls_dc, Rs_dc ] = get_CTMC_decorr_boundaries(sigmas, C, T, num_devs, sigma_dc);
Y_0s = [0 0 0];
    
% Form CTMC 1
center = Y_0s(1);
mu_func = @(s) drift_dc(1)*[s>-100000];
sig_func = @(s) sigma_dc(1)*[s>-100000];
[Q, y_1, c_index_1] = Q_Matrix(m_0, mu_func,sig_func,Ls_dc(1),Rs_dc(1),gridMethod,center, GridMultParam);
P1 = expm(Q*dt);

% Form CTMC 2
center = Y_0s(2);
mu_func = @(s) drift_dc(2)*[s>-100000];
sig_func = @(s) sigma_dc(2)*[s>-100000];
[Q, y_2, c_index_2] = Q_Matrix(m_0, mu_func,sig_func,Ls_dc(2),Rs_dc(2),gridMethod,center, GridMultParam);
P2 = expm(Q*dt);

% Form CTMC 3
center = Y_0s(3);
mu_func = @(s) drift_dc(3)*[s>-100000];
sig_func = @(s) sigma_dc(3)*[s>-100000];
[Q, y_3, c_index_3] = Q_Matrix(m_0, mu_func,sig_func,Ls_dc(3),Rs_dc(3),gridMethod,center, GridMultParam);
P3 = expm(Q*dt);


G = get_payoff_G_matrix_from_ygrid_3d( y_1, y_2, y_3, S_0s, sigmas, R, contractParams);

if contract == 1  % European
    % vals = exp(-r*T)*P1*G*P2.';
    vals = 0;
    for i=1:m_0
        for j=1:m_0
            for k=1:m_0
                vals = vals + P1(c_index_1,i)*P2(c_index_2,j)*P3(c_index_3,k)*G(i,j,k);
            end
        end
    end
    vals = vals*exp(-r*T);
    
end

end


