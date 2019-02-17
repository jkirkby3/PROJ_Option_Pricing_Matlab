function MNorm = MixedNormalRnd( n, p, a1,b1,a2,b2)
% Returns a column vector of size n of MixedNormal Random Variables
%   Detailed explanation goes here

Bern = rand(n,1) < p; %Column Vector of Bernoulli(p) RVs
Norm = randn(n,1);

MNorm = zeros(n,1);
MNorm(Bern ==1) = a1 +b1*Norm(Bern ==1);
MNorm(Bern ==0) = a2 +b2*Norm(Bern ==0);  

end

