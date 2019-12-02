function  price = PROJ_Swing_LinearRec( r,S_0,Dmax,rho_tau,T_0,T,Mtau,N,alpha,rnSYMB,Ks)

K1 = Ks(1); K2 = Ks(2); K3 = Ks(3); K4 = Ks(4);
w  = log(Ks/S_0);

Gx = @(x)G_func_swing( x,K1,K2,K3,K4,S_0 );

Dset = 1:Dmax; 
tau_RD = rho_tau*Dset;
tau1 = rho_tau; 

%%%%%%%%%%
K    = N/2;
dt   = tau1/Mtau;
nrdt = -r*dt;  

p = floor((T-T_0)/tau1);
Ttil_p = T - p*tau1;
Mtau_pr = floor(( Ttil_p - T_0)/dt);
M = p*Mtau + Mtau_pr;

xmin  = -alpha/2 + (w(3) + w(2))/2; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w      = log([K1,K2,K3,K4]./S_0);
dxtil  = 2*alpha/(N-1);  

nbars  = floor((w - xmin)/dxtil +1);
xnbars = xmin + dxtil*(nbars -1);

diffs = w - xnbars;
nbars(diffs<diffs(1)) = nbars(diffs<diffs(1)) - 1;

dx   = (w(4) - w(1))/(nbars(4) - nbars(1)); a = 1/dx;
xmin = w(1) - (nbars(1) - 1)*dx;

nbars(2:3)  = floor((w(2:3) - xmin)/dx +1);
xnbars      = xmin + dx*(nbars -1);
xnbars(1)   = w(1); xnbars(4) = w(4);

rhos   = w - xnbars;
zetastar  = a*rhos;
nnot   = floor(1-xmin*a);

xnot = xmin +(nnot-1)*dx;
xs = [xnot-dx,xnot,xnot+dx,xnot+2*dx];
inds = [nnot-1,nnot,nnot+1,nnot+2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  PHASE I %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

E    = S_0*exp(xmin + dx*(0:K-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  T^dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a2     = a^2;         zmin   = (1 - K)*dx;  %Kbar corresponds to zero
dw     = 2*pi*a/N;    DW     = dw*(1:N-1);
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

Cons4  = 1/12;  
xbars  = zeros(1,2);

Gs      = repmat(Gx(xmin + dx*(0:K-1))',1,Dmax);
THETG   = repmat(ThetaG,1,Dmax);

for d = 2:Dmax
    Gs(:,d) = d*Gs(:,d);
    THETG(:,d) = d*THETG(:,d);
end
G       = Gs(:,end);   %Dmax*G(x)
%%%----------------------------------------

edn    = exp(-dx);
dK21   = (K2 - K1);
dK43   = (K4 - K3);

Thetbar_dt = dK43*cumsum(beta(2*K:-1:K +1))' + dK21*[ fliplr(cumsum(beta(1:1:K-1)))';0];

for m = M-1:-1:M-(Mtau -1)
    pp      = ifft(toepM.*fft([THET(:,m+1);zeros(K,1)]));
    Cont_dt = pp(1:K) + Thetbar_dt;
    
    %%%----------------------------------------
    while Cont_dt(nms(1))> G(nms(1))
        nms(1) = nms(1) -1;
    end
    
    while Cont_dt(nms(2))> G(nms(2))
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

toepMDs = zeros(N,Dmax);

for d = 1:Dmax
    
    Cons3  = Cons1*exp(-r*tau_RD(d));
    grand  = grand1.*exp(tau_RD(d)*chfpoints);
    beta   = Cons3*real(fft([1/(24*a2) grand]));  
    toepMDs(:,d) = fft([beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)']);   %toepMDs(:,d) = fft(toepMDs(:,d));

end

Cont_Ds = zeros(K,Dmax);  
PSIs    = zeros(K,Dmax+1);  %Store Cont_dt in PSIs

ntil    = zeros(1,10);  %%% Just intialize at some number, say 10
dds     = zeros(1,10);


for m = M-Mtau:-1:1%m = M-Mtau :-1:1
    
    dstr = Dmax;
    while (m + dstr*Mtau) > M   %%% Find highest index of Cont_Ds which isn't zeros
        dstr = dstr - 1;
    end
   
    for d = 1:dstr
        pp      = ifft(toepMDs(:,d).*fft([THET(:,m+d*Mtau);zeros(K,1)]));  %make more efficient, dont keep recomputing m+d*Mtau
        Cont_Ds(:,d)  = pp(1:K) ;
        PSIs(:,d)     = Gs(:,d) + Cont_Ds(:,d);
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Now need to store Cont_dt in Dmax +1 spot
    pp      = ifft(toepM.*fft([THET(:,m+1);zeros(K,1)]));
    PSIs(:,Dmax +1) = pp(1:K) ;   %%%%Cont_dt

    %%% we only calculate PSIs up to dstr+1, but we still need to allow for exercise even when
	%%% the continuation value will be zero (also need to fix at end of algorithm)

    for d = dstr+1:Dmax
	PSIs(:,d) = Gs(:,d);
    end
    
    %%%PSIstr is PSI_Star
    [PSIstr,I]       = max(PSIs(:,1:Dmax+1),[],2);  %%% Note, since we only go up to dstr, returns index of dstr +1 instead of Dmax +1 when Cont_dt is maximizer
   
    I(nbars(2):nbars(3)) = Dmax+1;  %we know this holds in theory, ie continuation is optimal in the region of zero payoff
    
    [MstrL,ntil_L]   = max(PSIstr(1:nbars(1)));
    I(1:ntil_L)      = I(ntil_L);
    
    [MstrR,ntil_R]   = max(PSIstr(nbars(4):K));
    ntil_R           = nbars(4) - 1 + ntil_R;
    
    ntil_R = min(ntil_R,K-2);
    
    I(ntil_R:K)      = I(ntil_R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Find ntilJ
    
    ntil(1) = ntil_L;    dds(1)  = I(ntil(1) +1);
    THET(1:ntil(1) + 1,m)   = MstrL;   %% include the ntil_L+1 term for simplicity
    j = 2;
    while ntil(j-1) < ntil_R
        
        temp = find( I(ntil(j-1)+1 : ntil_R)~= dds(j-1),1 ) ;
        if isempty(temp)
            ntil(j) = ntil_R;
        else
            ntil(j) = min(ntil_R,temp + ntil(j-1) -1);
        end
        dds(j) = I(ntil(j) +1);
        
        x_n_til = xmin + dx*(ntil(j) -1);   %% ONLY Stores this one
        
        tmp1    = PSIs(ntil(j),dds(j-1)) - PSIs(ntil(j),dds(j));   
        tmp2    = PSIs(ntil(j)+1,dds(j-1)) - PSIs(ntil(j)+1,dds(j));
        
        if tmp1 ~= tmp2
            xtil = ((x_n_til +dx)*tmp1 - x_n_til*tmp2)/(tmp1-tmp2);  %%% we dont calculate these for ntil(1)
        else
            xtil = x_n_til;  
        end
        rho       = xtil - x_n_til;   
        zeta      = a*rho;  
        

        if dds(j-1) <=dstr %% PSI
            THET( ntil(j-1)+2 : ntil(j)-1,m) = THETG(ntil(j-1)+2 : ntil(j)-1,dds(j-1))...
                    + Cons4*( Cont_Ds(ntil(j-1)+1 : ntil(j)-2,dds(j-1)) + 10*Cont_Ds(ntil(j-1)+2 : ntil(j)-1,dds(j-1)) + Cont_Ds(ntil(j-1)+3 : ntil(j),dds(j-1)) );
        elseif dds(j-1) < Dmax +1      
            THET( ntil(j-1)+2 : ntil(j)-1,m) = THETG(ntil(j-1)+2 : ntil(j)-1,Dmax); %%% We know the best is Dmax
        elseif dds(j-1) == Dmax+1
            THET(ntil(j-1)+2 : ntil(j)-1,m) = Cons4*(PSIs(ntil(j-1)+1 : ntil(j)-2,dds(j-1)) + 10*PSIs(ntil(j-1)+2 : ntil(j)-1,dds(j-1)) + PSIs(ntil(j-1)+3 : ntil(j),dds(j-1)));
        end

        varths_dj         = Get_VarthsPSI_swing(zeta,ntil(j),PSIs(:,dds(j)),PSIs(:,dds(j-1)));
        THET(ntil(j),m)   = varths_dj(3) + varths_dj(1);
        THET(ntil(j)+1,m) = varths_dj(4) + varths_dj(2);   
        
        j = j+1;
    end
    J = j-1;
    
    THET(ntil_R +2: K, m) = MstrR;   %%% may need to change

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NEED TO VALUE AT END OF ALGORITHM (couldn't save THET(:,0))

dstr = Dmax;
while (dstr*Mtau) > M   %%% Find highest index of Cont_Ds which isn't zeros
    dstr = dstr - 1;
end

for d = 1:dstr
    pp      = ifft(toepMDs(:,d).*fft([THET(:,d*Mtau);zeros(K,1)]));  %make more efficient, dont keep recomputing m+d*Mtau
    Cont_Ds(:,d)  = pp(1:K) ;
    PSIs(:,d)     = Gs(:,d) + Cont_Ds(:,d);
end


pp      = ifft(toepM.*fft([THET(:,1);zeros(K,1)]));
PSIs(:,Dmax +1) = pp(1:K) ;   %%%%Cont_dt

for d = dstr+1:Dmax
    PSIs(:,d) = Gs(:,d);
end

[PSIstr,I]       = max(PSIs(:,1:Dmax+1),[],2);  %%% Note, since we only go up to dstr, returns index of dstr +1 instead of Dmax +1 when Cont_dt is maximizer

ys = PSIstr(inds);
price = spline(xs,ys,0);

end