function [ drift_dc, sigma_dc ] = decorrelate(sigmas, drifts, C, D )
% Applies the decorrelation transform
% C = C matrix, D = Diagonal Matrix
n = length(sigmas);
drift_dc = zeros(n,1);  
sigma_dc = zeros(n,1);

for i=1:n
    sigma_dc(i) = sqrt(D(i,i));
    sum_terms = 0;
    for j=1:n
        s = sigmas(j);
        sum_terms = sum_terms + (drifts(j) - s^2/2) * C(i,j)/s;
    end
    drift_dc(i) = sum_terms;
end


end

