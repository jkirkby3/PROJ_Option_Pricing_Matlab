function varths = Get_VarthsDD_swing(zetas,nms,Cont_D)
%
% 

zetas2     = zetas.^2; zetas3 = zetas.*zetas2; zetas4 = zetas.*zetas3;

varths = zeros(1,4);
gamms  = zeros(1,3);




gamms(1)  = (1 - zetas4(1))/8 +zetas3(1)/3 - zetas2(1)/4;
gamms(2)  = 5/12 + zetas4(1)/4 - zetas3(1)/3 - zetas2(1)/2 + zetas(1);
gamms(3)  = 1/12 - (1 + zetas4(1))/8 + zetas2(1)/4;

varths(3) = Cont_D(nms(1)-1)*gamms(1) + Cont_D(nms(1))*gamms(2) + Cont_D(nms(1)+1)*gamms(3);
%%%----------------------

gamms(1)  = zetas4(1)/8 - zetas3(1)/2 + zetas2(1)/2;
gamms(2)  = -zetas4(1)/4 + 2*zetas3(1)/3;
gamms(3)  = zetas4(1)/8 - zetas3(1)/6;

varths(4) = Cont_D(nms(1))*gamms(1) + Cont_D(nms(1)+1)*gamms(2) + Cont_D(nms(1)+2)*gamms(3);


%%%----------------------
%%% in this case varths_D(1) uses zetas(2)
gamms(1) = -1/24 +zetas4(2)/8 -zetas3(2)/3 +zetas2(2)/4;
gamms(2) = 5/12 - zetas4(2)/4 +zetas3(2)/3 +zetas2(2)/2 -zetas(2);
gamms(3) = (1 + zetas4(2))/8 - zetas2(2)/4;

varths(1) = Cont_D(nms(2)-1)*gamms(1) + Cont_D(nms(2))*gamms(2) + Cont_D(nms(2)+1)*gamms(3);
%%%----------------------
%%% in this case varths_D(2) uses zetas(2)
gamms(1)  = 1/12 - zetas4(2)/8 +(zetas3(2) - zetas2(2))/2;
gamms(2)  = 5/6  +zetas4(2)/4 - 2*zetas3(2)/3;
gamms(3)  = 1/12 -zetas4(2)/8 + zetas3(2)/6;

varths(2) = Cont_D(nms(2))*gamms(1) + Cont_D(nms(2)+1)*gamms(2) + Cont_D(nms(2)+2)*gamms(3);
%%%----------------------



end

