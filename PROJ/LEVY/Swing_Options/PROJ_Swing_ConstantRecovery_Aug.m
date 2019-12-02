function price = PROJ_Swing_ConstantRecovery_Aug( r,S_0,Dmax,T_0,T,tau1,Mtau,N,alpha,rnSYMB,Ks)
%  rnSYMB is risk-neutral Levy symbol
%  Ks = [K1,K2,K3,K4]

K1 = Ks(1); K2 = Ks(2); K3 = Ks(3); K4 = Ks(4);
w  = log(Ks/S_0);
xmin  = -alpha/2 + (w(3) + w(2))/2; 

Gx = @(x)G_func_swing( x,K1,K2,K3,K4,S_0 );

K    = N/2;
dt   = tau1/Mtau;
nrdt = -r*dt; 

p = floor((T-T_0)/tau1);
Ttil_p = T - p*tau1;
Mtau_pr = floor(( Ttil_p - T_0)/dt);

%Ttil_pp1 = Ttil_p -Mtau_pr*dt;  %If a final European step is needed
M = p*Mtau + Mtau_pr;

dxtil  = 2*alpha/(N-1);  

%================  Determine xGrid ========================================
nbars  = floor((w - xmin)/dxtil +1);
xnbars = xmin + dxtil*(nbars -1);

diffs = w - xnbars;
nbars(diffs<diffs(1)) = nbars(diffs<diffs(1)) - 1;

dx   = (w(4) - w(1))/(nbars(4) - nbars(1)); a = 1/dx;
xmin = w(1) - (nbars(1) - 1)*dx;

nbars(2:3)  = floor((w(2:3) - xmin)/dx +1);

rhos   = w - xnbars;
zetastar  = a*rhos;
nnot   = floor(1-xmin*a);
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Gaussian 3-point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;

%%%% PAYOFF CONSTANTS-----------------------------------
varthet_01   = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
zetas2     = zetastar.^2;
edn    = exp(-dx);

%%%----------------------------------------
%%% Initialize the dstars used in psis function
rhos_plus  = rhos*q_plus; rhos_minus = rhos*q_minus;
zetas_plus = a*rhos_plus; zetas_minus = a*rhos_minus;
eds1       = exp(rhos_minus); eds2 = exp(rhos/2); eds3 = exp(rhos_plus);

dbars_1    = zetas2/2;
dbars_0    = zetastar - dbars_1;         
ds_0       = zetastar.*(5*( (1-zetas_minus).*eds1 + (1-zetas_plus).*eds3 ) + 4*(2-zetastar).*eds2)/18;
ds_1       = edn*zetastar.*( 5*(zetas_minus.*eds1 + zetas_plus.*eds3) + 4*zetastar.*eds2 )/18; 

dstars    = zeros(1,4);
dstars(1) = dbars_0(2); dstars(2) = ds_0(2);
dstars(3) = ds_1(3); dstars(4) = dbars_1(3);

%==========================================================================
ThetaG    = GetThetaG_swing(xmin,K,dx,K1,K2,K3,K4,S_0 );

THET      = zeros(K,M);
THET(:,M) = Dmax*ThetaG;
E         = S_0*exp(xmin + dx*(0:K-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  PHASE I %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  T^dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a2     = a^2;
zmin   = (1 - K)*dx;  %Kbar corresponds to zero

dw     = 2*pi*a/N;
DW     = dw*(1:N-1);
grand1 = exp(-1i*zmin*DW ).*(sin(DW/(2*a))./DW ).^2./(2+cos(DW /a));
Cons1  = 24*a2/N;
%%%------------------------------------------------------------------
chfpoints = rnSYMB(DW);  

Cons2  = Cons1*exp(nrdt);
grand  = grand1.*exp(dt*chfpoints);
beta   = Cons2*real(fft([1/(24*a2) grand]));   %%%%  NOTE: all toep matrices incorporate exp(-r*dt)

toepM  = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];   toepM = fft(toepM);
%%%------------------------------------------------------------------

%%% initialize for search
nms    = zeros(1,2);
nms(1) = nbars(2) + 1;
nms(2) = nbars(3) - 1;
%%%----------------------------------------
Cons4  = 1/12;  
%%%----------------------------------------

xbars  = zeros(1,2);

%%%----------------------------------------
G      = Dmax*Gx(xmin + dx*(0:K-1))';
%%%----------------------------------------

edn    = exp(-dx);
dK21   = (K2 - K1);
dK43   = (K4 - K3);

Thetbar_dt = dK43*cumsum(beta(2*K:-1:K +1))' + dK21*[ fliplr(cumsum(beta(1:1:K-1)))';0];

for m =M-1:-1:M - (Mtau -1)
    pp      = ifft(toepM.*fft([THET(1:K,m+1);zeros(K,1)]));
    Cont_dt = pp(1:K) + Thetbar_dt;
   
    %%%----------------------------------------
    while (nms(1)>2) && (Cont_dt(nms(1))> G(nms(1)))
        nms(1) = nms(1) -1;
    end
    
    nms(2) = nms(2)+1;
    while nms(2) < K-2 && Cont_dt(nms(2))> G(nms(2))
        nms (2) = nms(2) + 1;
    end
    nms(2) = nms(2) -1;
    %%%----------------------------------------
    xnbars = xmin + dx*(nms -1);  
    
    %%%----------------------------------------
    tmp1     = Cont_dt(nms(1)) - G(nms(1));   tmp2 = Cont_dt(nms(1)+1) - G(nms(1)+1);
    xbars(1) = ((xnbars(1)+dx)*tmp1 - xnbars(1)*tmp2)/(tmp1-tmp2);
    %%%------------------

    tmp1     = Cont_dt(nms(2))-G(nms(2));     tmp2 = Cont_dt(nms(2)+1) - G(nms(2)+1);
    xbars(2) = ((xnbars(2)+dx)*tmp1 - xnbars(2)*tmp2)/(tmp1-tmp2);
    %%%------------------
    rhos       = xbars - xnbars;   
    zetas      = a*rhos;  
 
    psis      = Get_psis_swing_VER2( rhos,zetas, q_plus, q_minus, Ks, a , varthet_01, E,nms, nbars,edn,zetastar,dstars);
    varths_dt = Get_Varths_swing(zetas,nms,Cont_dt);

    %%%----------------------------------------
    THET(1:nms(1)-1,m) = THET(1:nms(1)-1,M);
    THET(nms(1),m)     = THET(nms(1),M) - Dmax*psis(1) +varths_dt(1);
    THET(nms(1)+1,m)   = Dmax*psis(2) + varths_dt(2);
    
    THET(nms(1)+2: nms(2)-1,m) = Cons4*(Cont_dt(nms(1)+1: nms(2)-2) + 10*Cont_dt(nms(1)+2: nms(2)-1) + Cont_dt(nms(1)+3: nms(2)));
    
    THET(nms(2),m)     = Dmax*psis(3) + varths_dt(3);
    THET(nms(2)+1,m)   = THET(nms(2)+1,M) - Dmax*psis(4) + varths_dt(4);
    THET(nms(2)+2:K,m) = THET(nms(2)+2:K,M);
    %%%----------------------------------------
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  PHASE II %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cons3  = Cons1*exp(-r*tau1);
grand  = grand1.*exp(tau1*chfpoints);
beta   = Cons3*real(fft([1/(24*a2) grand]));   
toepMD = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];   toepMD = fft(toepMD);

evec_D1 = cumsum(beta(2*K:-1:K +1))'; 
evec_D2 = [ fliplr(cumsum(beta(1:1:K-1)))';0];

count = Mtau;  %%%%starts the counter, when to reset the search indices

for m = M-Mtau:-1:1
    pp      = ifft(toepM.*fft([THET(1:K,m+1);zeros(K,1)]));
    Cont_dt = pp(1:K); %%%%%%%%%%%%%%%%%%
    
    pp      = ifft(toepMD.*fft([THET(1:K,m+Mtau);zeros(K,1)]));
    Cont_D  = pp(1:K) + evec_D2*THET(1,m+Mtau) + evec_D1*THET(K,m+Mtau);
    PSI     = G + Cont_D;
    
    %%%% Reset at each of the discontinuities in the value function
    if rem(count,Mtau)==0
        nms(1) = nbars(2); nms(2) = nbars(3)+1;
    end
    count = count+1;
    %%%----------------------------------------
    while nms(1)>2 && Cont_dt(nms(1))> PSI(nms(1))
        nms(1) = nms(1) -1;
    end
    
    while nms(2)<K-2 && Cont_dt(nms(2))> PSI(nms(2))
        nms (2) = nms(2) + 1;

    end

     nms(2) = nms(2)-1;
    %%%----------------------------------------
    xnbars = xmin + dx*(nms -1);
    
    %%%----------------------------------------
    tmp1     = Cont_dt(nms(1))-PSI(nms(1));        tmp2 = Cont_dt(nms(1)+1) - PSI(nms(1)+1);
    xbars(1) = xnbars(1) + max(0,dx*tmp1/(tmp1-tmp2));
    
    %%%------------------
    tmp1     = Cont_dt(nms(2))-PSI(nms(2));        tmp2 = Cont_dt(nms(2)+1) - PSI(nms(2)+1);
    xbars(2) = xnbars(2) + max(0,dx*tmp1/(tmp1-tmp2));
      
    %%%------------------
    rhos      = xbars - xnbars;   
    zetas     = a*rhos ;
    %%%----------------------------------------
    psis        = Get_psis_swing_VER2( rhos,zetas, q_plus, q_minus, Ks, a , varthet_01, E,nms, nbars,edn,zetastar,dstars);
    
    varths_dt  = Get_Varths_swing(zetas,nms,Cont_dt);
    varths_D   = Get_VarthsDD_swing(zetas,nms,Cont_D);
    
    
    %%%----------------------------------------
    THET(2:nms(1)-1,m) = THET(2:nms(1)-1,M) + Cons4*(Cont_D(1:nms(1)-2) + 10*Cont_D(2:nms(1)-1) + Cont_D(3:nms(1)));
    THET(1,m)          = THET(2,m);   %%%% For simplicity... won't affect solution
    
    THET(nms(1),m)     = THET(nms(1),M) - Dmax*psis(1) + varths_dt(1)    + varths_D(3);
    THET(nms(1)+1,m)   = Dmax*psis(2) + varths_dt(2)    + varths_D(4);
    
    THET(nms(1)+2: nms(2)-1,m) = Cons4*(Cont_dt(nms(1)+1: nms(2)-2) + 10*Cont_dt(nms(1)+2: nms(2)-1) + Cont_dt(nms(1)+3: nms(2)));
    
    THET(nms(2),m)       = Dmax*psis(3) + varths_dt(3)     + varths_D(1);
    THET(nms(2)+1,m)     = THET(nms(2)+1,M) - Dmax*psis(4) + varths_dt(4)     + varths_D(2);
    
    THET(nms(2)+2:K-1,m) = THET(nms(2)+2:K-1,M)  + Cons4*(Cont_D(nms(2)+1:K-2) + 10*Cont_D(nms(2)+2:K-1) + Cont_D(nms(2)+3:K));
    THET(K,m)            = THET(K-1,m); %%%% For simplicity... won't affect solution
    %%%----------------------------------------
end

pp      = ifft(toepM.*fft([THET(1:K,1);zeros(K,1)]));
Cont_dt = pp(1:K);

pp      = ifft(toepMD.*fft([THET(1:K,Mtau);zeros(K,1)]));
Cont_D  = pp(1:K) + evec_D2*THET(1,Mtau) + evec_D1*THET(K,Mtau);
PSI     = G + Cont_D;

xnot = xmin +(nnot-1)*dx;

xs   = [xnot-dx,xnot,xnot+dx,xnot+2*dx];
inds = [nnot-1,nnot,nnot+1,nnot+2];
ys   = max(PSI(inds),Cont_dt(inds));

price = spline(xs,ys,0);

end

