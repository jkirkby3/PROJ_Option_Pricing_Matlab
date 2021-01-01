function price = PROJ_Swing_FixedRights(M, N, alpha, rnCHF, r, Dmax, T_0, T, Ns, S_0, Ks)
%  
%   FOR NOW: T_0 must = 0
%

K = N/2;
w  = log(Ks/S_0);
K1 = Ks(1); K2 = Ks(2); K3 = Ks(3); K4 = Ks(4);

Gx = @(x)G_func_swing( x,K1,K2,K3,K4,S_0 );

N    = 2*K; 
dt = (T-T_0)/M;
nrdt = -r*dt; 
xmin  = -alpha/2 + (w(3) + w(2))/2; 

dxtil  = 2*alpha/(N-1);  

nbars  = floor((w - xmin)/dxtil +1);
xnbars = xmin + dxtil*(nbars -1);

diffs = w - xnbars;
nbars(diffs<diffs(1)) = nbars(diffs<diffs(1)) - 1;

dx   = (w(4) - w(1))/(nbars(4) - nbars(1)); a = 1/dx;
xmin = w(1) - (nbars(1) - 1)*dx;

nbars(2:3)  = floor((w(2:3) - xmin)/dx +1);
nnot   = floor(1-xmin*a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  PHASE I %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%----------------------------------------
G      = Dmax*Gx(xmin + dx*(0:K-1))';
%%%----------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Gaussian 3-point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;

%%%% PAYOFF CONSTANTS-----------------------------------
varthet_01   = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;

ThetaDmax    = Dmax*GetThetaG_swing(xmin,K,dx,K1,K2,K3,K4,S_0 );
THET         = repmat(ThetaDmax,1,Ns+1);   %%% For fixed rights case, dont need to store for differnt time levels
E            = S_0*exp(xmin + dx*(0:K-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Cons2  = Cons1*exp(nrdt);
grand  = grand1.*rnCHF(DW );
beta   = Cons2*real(fft([1/(24*a2) grand]));   %%%%  NOTE: all toep matrices incorporate exp(-r*dt)
toepM  = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];   toepM = fft(toepM);
%%%------------------------------------------------------------------

%%% initialize for search
nms    = zeros(1,2);
nms(1) = nbars(2) + 1;
nms(2) = nbars(3) ;
%%%----------------------------------------
Cons4  = 1/12;  
%%%----------------------------------------
xbars  = zeros(1,2);
edn    = exp(-dx);

nm1vec = ones(Ns+1,1)*nbars(2);
nm2vec = ones(Ns+1,1)*nbars(3);

CONT_nm1 = zeros(K,1);

for m =M-1:-1:1
    
    for n = 2:Ns+1   %  the first column is just zeros
        pp           = ifft(toepM.*fft([THET(1:K,n);zeros(K,1)]));
        CONT_n     = pp(1:K);
        PSI          = G + CONT_nm1;
        
        %%%----------------------------------------
        if n-1 > M -m %in these states it is always optimal to exercise
            THET(:,n) = ThetaDmax;
            THET(2:K-1,n) = THET(2:K-1,n) + Cons4*(CONT_nm1(1:K-2) + 10*CONT_nm1(2:K-1) + CONT_nm1(3:K));
            nm1vec(n) = nbars(2); nm2vec(n) = nbars(3);
        else %n -1 <= M - m  %we added 1 to each state n
            
             nms(1) = nbars(2); nms(2) = nbars(3) +1; %%% NOTE: this is more accurate for very low resolutions
            %nms(1) = nm1vec(n); nms(2) = nm2vec(n) +1;  %%%NOTE: this is slower than reseting based on previous EE point
            
            while nms(1)>1 && CONT_n(nms(1))> PSI(nms(1))  %%% If they are equal to start, then nms(1) = nbars(2)
                nms(1) = nms(1) -1;
            end

            while nms(2)<K-1 && CONT_n(nms(2))> PSI(nms(2)) %%%% less than K-1 so the algorithm works in extreme cases
                nms(2) = nms(2) + 1;
            end
            nms(2) = nms(2) -1;
            
            nm1vec(n) = nms(1); nm2vec(n) = nms(2);

            xnbars = xmin + dx*(nms -1);

            %%%----------------------------------------
            tmp1     = CONT_n(nms(1))-PSI(nms(1));      tmp2 = CONT_n(nms(1)+1) - PSI(nms(1)+1);
            xbars(1) = xnbars(1) +max(0,dx*tmp1/(tmp1-tmp2));
            
            
            %%%------------------
            tmp1     = CONT_n(nms(2))-PSI(nms(2));      tmp2 = CONT_n(nms(2)+1) - PSI(nms(2)+1);
            xbars(2) = xnbars(2) +max(0,dx*tmp1/(tmp1-tmp2));
        
            %%%------------------
            rhos      = xbars - xnbars;
            zetas     = a*rhos;
            %%%----------------------------------------

            psis          = Get_psis_swing( rhos,zetas, q_plus, q_minus, Ks, a , varthet_01, E,nms, nbars,edn);
            varths_dt_n   = Get_Varths_swing(zetas,nms,CONT_n);
            varths_dt_nm1 = Get_VarthsDD_swing(zetas,nms,CONT_nm1);
            
            %%%----------------------------------------
            THET(1,n)          = PSI(1);   %%%% For simplicity... won't affect solution
            THET(2:nms(1)-1,n) = ThetaDmax(2:nms(1)-1) + Cons4*(CONT_nm1(1:nms(1)-2) + 10*CONT_nm1(2:nms(1)-1) + CONT_nm1(3:nms(1)));

            THET(nms(1),n)     = ThetaDmax(nms(1)) - Dmax*psis(1) + varths_dt_n(1)    + varths_dt_nm1(3);
            THET(nms(1)+1,n)   = Dmax*psis(2) + varths_dt_n(2)    + varths_dt_nm1(4);

            THET(nms(1)+2: nms(2)-1,n) = Cons4*(CONT_n(nms(1)+1: nms(2)-2) + 10*CONT_n(nms(1)+2: nms(2)-1) + CONT_n(nms(1)+3: nms(2)));

            THET(nms(2),n)       = Dmax*psis(3) + varths_dt_n(3)     + varths_dt_nm1(1);
            THET(nms(2)+1,n)     = ThetaDmax(nms(2)+1) - Dmax*psis(4) + varths_dt_n(4)     + varths_dt_nm1(2);

            THET(nms(2)+2:K-1,n) = ThetaDmax(nms(2)+2:K-1)  + Cons4*(CONT_nm1(nms(2)+1:K-2) + 10*CONT_nm1(nms(2)+2:K-1) + CONT_nm1(nms(2)+3:K));
            THET(K,n)            = PSI(K); %%%% For simplicity... won't affect solution

        end
        
        CONT_nm1 = CONT_n;
    end
   CONT_nm1 = zeros(K,1);
end

%%%%%%%%%%   Time Zero  %%%%%%%%%%%%
%%% Have N_s rights
pp           = ifft(toepM.*fft([THET(1:K,Ns);zeros(K,1)]));
CONT_nm1     = pp(1:K);

pp           = ifft(toepM.*fft([THET(1:K,Ns+1);zeros(K,1)]));
CONT_n     = pp(1:K);


PSI = G + CONT_nm1;  %%% correspons to Ns-1

xnot = xmin +(nnot-1)*dx;
xs   = [xnot-dx,xnot,xnot+dx,xnot+2*dx];
inds = [nnot-1,nnot,nnot+1,nnot+2];
ys   = max(PSI(inds),CONT_n(inds));

price = spline(xs,ys,0);

end

