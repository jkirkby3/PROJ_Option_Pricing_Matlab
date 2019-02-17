function Poi = PoissonRnd( n, theta)
% Returns a column vector of size n of Poisson Random Variables
% For a Poisson process, use theta = lambda*dt, where lambda is arrival rate
% Density: exp(-theta)*(theta)^k/(k!)

%Uses method on page 128, Figure 3.9 of glasserman
Unif = rand(n,1);
Poi = zeros(n,1);

for j = 1:n
    p = exp(-theta); 
    F = p;
    N = 0;
    U = Unif(j);
    while U > F
        N = N + 1;
        p = p*theta/N;
        F = F + p;
    end
    Poi(j) = N;
end

end

