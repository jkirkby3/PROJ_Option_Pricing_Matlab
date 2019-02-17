function DExpo = DoubleExpoRnd(n, p_up, eta1,eta2)
% Returns a column vector of size n of DoubleExponetial Random Variables
% Note: To Generate an expo(mean=1/lambda), use E = -log(Unif)/lambda
% eta1 is expo param for up jumps, eta2 for down jumps

Bern = rand(n,1) < p_up; %Column Vector of Bernoulli(p_up) RVs

Unif  = rand(n,1);
DExpo = zeros(n,1);
DExpo(Bern ==1) = -log(Unif(Bern ==1))/eta1;
DExpo(Bern ==0) = log(Unif(Bern ==0))/eta2;  %negative jumps, -E

end

