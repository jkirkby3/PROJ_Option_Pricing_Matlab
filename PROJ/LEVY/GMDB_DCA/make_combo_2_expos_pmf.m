function p = make_combo_2_expos_pmf( b1, b2, xi1, xi2, Nmax )
% Creates a pmf from combination of 2 exponentials

p = b1*exp(-xi1*(1:Nmax))*(exp(xi1) - 1) + b2*exp(-xi2*(1:Nmax))*(exp(xi2) - 1);

temp = b1*(1 - exp(-xi1*(Nmax - 1))) + b2*(1 - exp(-xi2*(Nmax - 1)));
p(end) = 1 - temp;

end

