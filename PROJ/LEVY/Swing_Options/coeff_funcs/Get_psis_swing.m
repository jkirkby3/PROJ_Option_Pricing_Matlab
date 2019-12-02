function psis = Get_psis_swing( rhos,zetas, q_plus, q_minus, Ks, a , varthet_01, E,nms, nbars,edn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
K1 = Ks(1); K2 = Ks(2); K3 = Ks(3); K4 = Ks(4);
dK21 = K2 - K1; dK43 = K4- K3;

psis = zeros(1,4);

zetas2     = zetas.^2;

%%%----------------------------------------
rhos_plus  = rhos*q_plus; rhos_minus = rhos*q_minus;
zetas_plus = a*rhos_plus; zetas_minus = a*rhos_minus;
eds1       = exp(rhos_minus); eds2 = exp(rhos/2); eds3 = exp(rhos_plus);
dbars_1    = zetas2/2;
dbars_0    = zetas - dbars_1;         %dbars_1 = zetas + .5*((zetas - 1).^2 - 1);
ds_0       = zetas.*(5*( (1-zetas_minus).*eds1 + (1-zetas_plus).*eds3 ) + 4*(2-zetas).*eds2)/18;
ds_1       = edn*zetas.*( 5*(zetas_minus.*eds1 + zetas_plus.*eds3) + 4*zetas.*eds2 )/18; 
%%%-------------------

if nms(1) >= nbars(1)
    if nms(1)<nbars(2)
        psis(1) = K2*(.5 - dbars_0(1)) - E(nms(1))*(varthet_01 - ds_0(1));
    else %if nms(1)=nbars(2)
        psis(1) = 0;
    end
    psis(2) = K2*dbars_1(1) - E(nms(1)+1)*ds_1(1);
else
    psis(1) = dK21*(.5 - dbars_0(1));   %dK21 = (K2 - K1)
    psis(2) = dK21*dbars_1(1);
end

%%%-------------------
if nms(2) < nbars(4)
    psis(3) = E(nms(2))*(varthet_01 - ds_0(2)) - K3*(.5 - dbars_0(2));
    if nms(2)>nbars(3)
        psis(4) = E(nms(2)+1)*ds_1(2) - K3*dbars_1(2);
    else
        psis(4) = 0;
    end
else
    psis(3) = dK43*(.5 - dbars_0(2));   %dK43 = (K4 - K3)
    psis(4) = dK43*dbars_1(2);
end
%%%----------------------------------------


end

