function y = G_func_swing( x,K1,K2,K3,K4,S_0 )
% K1,..,K4 are the kink points
% G(x) is a function of x = ln(S/S_0)
w1 = log(K1/S_0); w2 = log(K2/S_0); w3 = log(K3/S_0); w4 = log(K4/S_0);

y = zeros(size(x));

y(x <= w1) = K2 - K1;
y(w1 < x & x <= w2) =  K2 - S_0*exp(x(w1 < x & x <= w2));
y(x < w4 & x > w3) = S_0*exp(x(x < w4 & x > w3)) - K3;
y(x >= w4) = K4 - K3;

   
end

