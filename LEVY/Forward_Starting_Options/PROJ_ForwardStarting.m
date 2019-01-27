function price = PROJ_ForwardStarting(N,alph,r,q,T1,T2,S_0,call,rnCHF1,rnCHF2)
% For now Order must = 3
%   Detailed explanation goes here
T = T1 + T2;

dx = 2*alph/(N-1); a = 1/dx;
dw = 2*pi/(N*dx);

omega = dw*(1:N-1);  %We calcuate coefficient of w=0 explicitly
nbar = N/2;
xmin = (1 - N/2)*dx;

%%%% Cubic Spline
b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
grand = @(w)rnCHF2(w).*(sin(w/(2*a))./w).^4./(b0 + b1*cos(w/a) +b2*cos(2*w/a) +b3*cos(3*w/a));
beta  = real(fft([1/(32*a^4) exp(-1i*xmin*omega).*feval(grand,omega)]));

%FIND Value of E[(1 - exp(X_tau))^+]
G = zeros(1,N); 

G(nbar +1) = 1*(1/24 - 1/20*exp(dx)*(exp(-7/4*dx)/54 + exp(-1.5*dx)/18 + exp(-1.25*dx)/2 + 7*exp(-dx)/27));

G(nbar )   =  1*(.5 -.05*(28/27 + exp(-7/4*dx)/54 + exp(-1.5*dx)/18 + exp(-1.25*dx)/2 + 14*exp(-dx)/27 ...
            + 121/54*exp(-.75*dx) + 23/18*exp(-.5*dx) + 235/54*exp(-.25*dx)));

G(nbar -1) = 1*( 23/24 - exp(-dx)/90*( (28 + 7*exp(-dx))/3 ...
            + ( 14*exp(dx) + exp(-7/4*dx) + 242*cosh(.75*dx) + 470*cosh(.25*dx))/12 ...
            +.25*(exp(-1.5*dx) + 9*exp(-1.25*dx) + 46*cosh(.5*dx))) );

vartheta_star = 1/90*( 14/3*(2+cosh(dx)) ...
                + .5*(cosh(1.5*dx) + 9*cosh(1.25*dx) +23*cosh(.5*dx))...
                +  1/6*(cosh(7/4*dx) + 121*cosh(.75*dx) +235*cosh(.25*dx)));        

G(1: nbar -2) = 1 - 1*exp(xmin +dx*(0:nbar-3))*vartheta_star;
Cons = 32*a^4;

Vbar2 = Cons/N*G(1,1:(nbar +1))*(beta(1,1:(nbar +1))');

%%% Find Second Expansion
grand = @(w)rnCHF1(w).*(sin(w/(2*a))./w).^4./(b0 + b1*cos(w/a) +b2*cos(2*w/a) +b3*cos(3*w/a));
beta  = real(fft([1/(32*a^4) exp(-1i*xmin*omega).*feval(grand,omega)]));

G(1:N) = exp(xmin + (0:(N-1))*dx);
Vbar1 = Cons/N*G*beta';
Val_put = exp(-r*T)*Vbar2*Vbar1*S_0*vartheta_star;


%%%%  PUT CALL PARITY/PRICING FORMULA
if call ==1
    Val_Proj = S_0*(exp(-q*T)-exp(-r*T2)*exp(-q*T1)) + Val_put;  %check this formula
else
    Val_Proj = Val_put;
end

price = Val_Proj;

end

