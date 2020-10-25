function [paths_1, paths_2] = Simulate_Diffusion_2D(S_0s, drifts, sigmas, rho,  N_sim, M, dt, exponential)
% R = correlation matrix
% 

paths_1 = zeros(N_sim,M+1);
paths_2 = zeros(N_sim,M+1);

paths_1(:,1) = S_0s(1);
paths_2(:,1) = S_0s(2);

if exponential == 1  % Need to convexity adjust the drift
    drifts = (drifts - .5*sigmas.^2);
end
drifts = drifts * dt;

Sigsqdt = sigmas*sqrt(dt);
rhosq = sqrt(1 - rho^2);

if exponential
    for m = 1:M
        W1 = randn(N_sim,1); 
        W2 = randn(N_sim,1); 

        paths_1(:,m+1) = paths_1(:,m).*exp(drifts(1) + Sigsqdt(1)*W1);  %log scheme
        paths_2(:,m+1) = paths_2(:,m).*exp(drifts(2) + Sigsqdt(2)*(rho*W1 + rhosq*W2));  %log scheme
    end
else
    for m = 1:M
        W1 = randn(N_sim,1); 
        W2 = randn(N_sim,1); 

        paths_1(:,m+1) = paths_1(:,m) + drifts(1) + Sigsqdt(1)*W1;  %log scheme
        paths_2(:,m+1) = paths_2(:,m) + drifts(2) + Sigsqdt(2)*(rho*W1 + rhosq*W2);  %log scheme
    end
end
    

end

