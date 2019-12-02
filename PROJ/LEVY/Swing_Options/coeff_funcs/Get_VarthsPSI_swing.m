function varths = Get_VarthsPSI_swing(zetaj,ntilj,PSIdj,PSIdjm1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


zetaj2     = zetaj^2; zetaj3 = zetaj*zetaj; zetaj4 = zetaj*zetaj3;

varths = zeros(1,4);
gamms  = zeros(1,3);


gamms(1)  = (1 - zetaj4)/8 + zetaj3/3 - zetaj2/4;
gamms(2)  = 5/12 + zetaj4/4 - zetaj3/3 - zetaj2/2 + zetaj;
gamms(3)  = 1/12 - (1 + zetaj4)/8 + zetaj2/4;

varths(3) = PSIdjm1(ntilj(1)-1)*gamms(1) + PSIdjm1(ntilj(1))*gamms(2) + PSIdjm1(ntilj(1)+1)*gamms(3);

gamms(1) = 1/12 - gamms(1);
gamms(2) = 5/6  - gamms(2);
gamms(3) = 1/12 - gamms(3);

varths(1) = PSIdj(ntilj(1)-1)*gamms(1) + PSIdj(ntilj(1))*gamms(2) + PSIdj(ntilj(1)+1)*gamms(3);

%%%----------------------

gamms(1)  = zetaj4/8 - zetaj3/2 + zetaj2/2;
gamms(2)  = -zetaj4/4 + 2*zetaj3/3;
gamms(3)  = zetaj4/8 - zetaj3/6;

varths(4) = PSIdjm1(ntilj(1))*gamms(1) + PSIdjm1(ntilj(1)+1)*gamms(2) + PSIdjm1(ntilj(1)+2)*gamms(3);

gamms(1) = 1/12 - gamms(1);
gamms(2) = 5/6  - gamms(2);
gamms(3) = 1/12 - gamms(3);

varths(2) = PSIdj(ntilj(1))*gamms(1) + PSIdj(ntilj(1)+1)*gamms(2) + PSIdj(ntilj(1)+2)*gamms(3);

end

