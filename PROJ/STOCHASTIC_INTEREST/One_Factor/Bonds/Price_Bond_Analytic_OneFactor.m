function [ price, r_mean, r_var] = Price_Bond_Analytic_OneFactor( T, model, modparam )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r0 = modparam.r0;

if model == 1
    
    k = modparam.eta; theta = modparam.theta; sig = modparam.Sigmar; 
    sig2 = sig^2;
    h = sqrt(k^2 + 2*sig2);
    
    num1 = 2*h*exp((k+h)*T/2);
    num2 = 2*h + (k+h)*(exp(T*h) - 1);
    A = (num1/num2)^(2*k*theta/sig2);
    
    num1 = 2*(exp(T*h) - 1);

    B = num1/num2;
    price = A*exp(-B*r0);  % Affine model
    
    ekT = exp(-k*T);
    r_mean = r0*ekT + theta*(1 - ekT);
    r_var = sig2/k*( r0*(ekT - exp(-2*k*T)) + theta/2 * (1 - ekT)^2);
    
elseif model == 2
    eta = modparam.eta; theta = modparam.theta; sig = modparam.Sigmar; k = eta;
    sig2 = sig^2;
    
    B = (1-exp(-k*T))/k;
    A = exp((theta - sig2/(2*k*k))*(B - T) - sig2/(4*k)*B*B);
    price = A*exp(-B*r0);  % Affine model
    
    r_mean = exp(-eta*T)*r0 + theta*(1-exp(-eta*T));
    r_var = sig2/(2*eta)*(1 -exp(-2*eta*T));
end

end

