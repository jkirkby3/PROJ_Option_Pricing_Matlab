function theta = GetThetaG_swing(xmin,K,dx,K1,K2,K3,K4,S_0 )
%
% 
a = 1/dx;
w      = log([K1,K2,K3,K4]./S_0);
nbars  = floor(a*(w - xmin) +1);
xnbars = xmin + dx*(nbars -1);
rhos   = (w - xnbars);
zetas  = a*rhos;

theta = zeros(K,1);
theta(1:nbars(1) - 1) = K2 - K1;
theta(nbars(4) +2:K) = K4 - K3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Gaussian 3-point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;

%%%% PAYOFF CONSTANTS-----------------------------------
varthet_01 = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
varthet_m10 = exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4)/18;
varthet_star = varthet_01 + varthet_m10;


E = S_0*exp(xmin + dx*(0:K-1));

zetas2     = zetas.^2; zetas3 = zetas.*zetas2; zetas4 = zetas.*zetas3;
rhos_plus  = rhos*q_plus; rhos_minus = rhos*q_minus;
zetas_plus = a*rhos_plus; zetas_minus = a*rhos_minus;
eds1 = exp(rhos_minus); eds2 = exp(rhos/2); eds3 = exp(rhos_plus);

dbars_1 = zetas2/2;
dbars_0 = zetas - dbars_1;         %dbars_1 = zetas + .5*((zetas - 1).^2 - 1);

ds_0    = zetas.*(5*( (1-zetas_minus).*eds1 + (1-zetas_plus).*eds3 ) + 4*(2-zetas).*eds2)/18;
ds_1    = exp(-dx)*zetas.*( 5*(zetas_minus.*eds1 + zetas_plus.*eds3) + 4*zetas.*eds2 )/18; 

%%% Get Intial Coeffs
theta(1:nbars(1) - 1) = K2 - K1;
theta(nbars(4) +2:K)  = K4 - K3;

theta(nbars(1)) = K2 - K1*(.5 + dbars_0(1)) - E(nbars(1))*(varthet_01 - ds_0(1));
theta(nbars(1)+1) = K2 - K1*dbars_1(1) - E(nbars(1)+1)*(varthet_star - ds_1(1));

theta(nbars(1)+2:nbars(2)-1) = K2 - varthet_star*E(nbars(1)+2:nbars(2)-1);

tmp1 = K2*(.5 + dbars_0(2) - exp(-rhos(2))*(varthet_m10 + ds_0(2)));
tmp2 = K2*(dbars_1(2) - exp(dx-rhos(2))*ds_1(2));

tmp3 = K3*( (varthet_01 - ds_0(3))*exp(-rhos(3)) - (.5 - dbars_0(3)));
tmp4 = K3*(exp(dx-rhos(3))*(varthet_star - ds_1(3)) - (1-dbars_1(3)));

if K3>K2
    theta(nbars(2))   = tmp1;
    theta(nbars(2)+1) = tmp2;
    theta(nbars(3))   = tmp3;
    theta(nbars(3)+1) = tmp4;
elseif K3==K2
    theta(nbars(2))   = tmp1 +tmp3;
    theta(nbars(2)+1) = tmp2 + tmp4;
end

theta(nbars(3) +2:nbars(4) -1) = E(nbars(3) +2:nbars(4) -1)*varthet_star - K3;

theta(nbars(4)) = E(nbars(4))*(varthet_m10 + ds_0(4)) - K3 + K4*(.5 -dbars_0(4));
theta(nbars(4)+1) = K4*(1 - dbars_1(4)) - K3 + E(nbars(4)+1)*ds_1(4);


end

