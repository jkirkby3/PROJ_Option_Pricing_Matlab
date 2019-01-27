function price = ParisianPROJ_alpha(N, call, down, S_0, W, H, M, r, rnCHF, T,Gamm, resetting, alph)
% S_0 = Initial price
% call = 1 for call (otherwise it's a put)
% down = 1 for down and out (otherwise it's up and out)
% N = number of basis elements (power of 2, e.g. 2^12)
% alph = size of valuation grid (width of projected density support is ~2*alph) .. fix based on density(T), rather than density(dt)
% H    = barrier
% rebate = rebate paid immediately upon passing the barrier (knocking-out)
% r = Interest rate
% q = dividend yield
% T = Time (in years)
% M = number of discrete monitoring points
% rnCHF = risk netural characteristic function
% Gamm = maximum number of discretely monitored excursions into the knockout region allowed (more than this results in knockout)
% resetting = 1 if a reseting type parisian option (otherwise its cumulative, ie never resets)

if ~(down == 1 && call ~= 1) && ~(down ~=1 && call == 1)
    fprintf('Sorry, currently only Up and out calls, and down and out puts have been coded \n')
    return
end

gamm0 = 1; %HARDCODED: this param would allow us to specify an inital consumed budget

dt   = T/M;
nrdt = -r*dt;
h = log(H/S_0);
lws   = log(W/S_0);

%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;

dx = 2*alph/(N-1); 
xmin = -alph/2;

n_h = floor((h-xmin)/dx +1); 
xmin = h - (n_h -1)*dx;    %realign so that h is on the grid (this is important for the case where h = 0)


if h~= 0 %Realign so that h and 0 are both members of grid (if possible)
    nnot =  floor(1-xmin/dx);
    if abs(h) > dx  %so that n_h ~= nnot
        dx = (h - 0)/(n_h - nnot);
        xmin = dx*(1-nnot);  %hence nnot should remain on the grid
        %n_h = floor((h-xmin)/dx +1);  %NOT Numerically Stable
        n_h = floor(nnot + h/dx);  %Numerically Stable
    end
else 
    nnot = n_h;  
end

a    = 1/dx;
a2   = a^2;
zmin = (1 - N/2)*dx;  %Kbar corresponds to zero

%Cons = 24*a2/N;
Cons2 = 24*a2*exp(nrdt)/N;
dw    = 2*pi*a/N;
grand = dw*(1:N-1);
grand = exp(-1i*zmin*grand).*rnCHF(grand).*(sin(grand/(2*a))./grand).^2./(2+cos(grand/a));
beta  = Cons2*real(fft([1/(24*a2) grand]));   %%%%  NOTE: all toep matrices incorporate exp(-r*dt)


interp_Atend = 0;
if 0 < abs(h) && abs(h)<dx
    interp_Atend = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   DETERMINE COMMON Params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K     = N/2;
nbar  = floor(a*(lws - xmin)+1);
rho   = lws - (xmin+(nbar - 1)*dx);
zeta  = a*rho;

toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];   toepM = fft(toepM);

%%%% PAYOFF CONSTANTS-----------------------------------
varthet_01 = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
varthet_m10 = exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4)/18;
varthet_star = varthet_01 + varthet_m10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if down == 1 && call ~= 1  %DOP
    %l = log(H/S_0);
    
    n_l = n_h;

    zeta_plus = zeta*q_plus; zeta_minus = zeta*q_minus;
    rho_plus = rho*q_plus; rho_minus = rho*q_minus;

    ed1 = exp(rho_minus); ed2 = exp(rho/2); ed3 = exp(rho_plus);

    dbar_1 = zeta^2/2;
    dbar_0 = zeta - dbar_1;         %  dbar_1 = zeta + .5*((zeta - 1)^2 - 1);
    d_0    = zeta*(5*( (1-zeta_minus)*ed1 + (1-zeta_plus)*ed3 ) + 4*(2-zeta)*ed2)/18;
    d_1    = zeta*( 5*(zeta_minus*ed1 + zeta_plus*ed3) + 4*zeta*ed2 )/18;            

    %%%%Thet(1)        =  W/2 - H*varthet_01;
    Thet = zeros(K,1);
    Thet(1:nbar-1) =  W - exp(xmin +dx*(0:nbar-2))*S_0*varthet_star;
    Thet(nbar)     =  W*(.5 + dbar_0 - exp(-rho)*(varthet_m10 + d_0));
    Thet(nbar + 1) =  W*(dbar_1 - exp(- rho)*d_1);
  
    %%%%%%% Initialize Val
    Val = zeros(K,Gamm+1);
    p   = ifft(toepM.*fft([Thet;zeros(K,1)]));
    for j=1:Gamm
       Val(:,j) = p(1:K); 
    end
    Thet(1:n_l-1) = 0;  %up to the barrier l, the coefficients are zero
    Thet(n_l) = W/2 - H*varthet_01;  %because this is the case where you can knock out at termination
    p   = ifft(toepM.*fft([Thet;zeros(K,1)]));
    Val(:,Gamm+1) = p(1:K);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% RESETTING PARISIAN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if resetting == 1
        for m=M-2:-1:0
            %These Thet must be kept outside of loop, else they will be inadvertently altered
            Thet(n_l+1:K-1) = (Val(n_l:K-2,1)+10*Val(n_l+1:K-1,1)+Val(n_l+2:K,1))/12;
            Thet(K)         = (13*Val(K,1)+15*Val(K-1,1)-5*Val(K-2,1)+Val(K-3,1))/48;
            ThetPartial     = (13*Val(n_l,1)+15*Val(n_l+1,1)-5*Val(n_l+2,1)+Val(n_l+3,1))/48;
            for j=1:Gamm

                Thet(1)      = (13*Val(1,j+1)+15*Val(2,j+1)-5*Val(3,j+1)+Val(4,j+1))/48;
                Thet(2:n_l-1) = (Val(1:n_l-2,j+1)+10*Val(2:n_l-1,j+1)+Val(3:n_l,j+1))/12;
                Thet(n_l)     = (13*Val(n_l,j+1)+15*Val(n_l-1,j+1)-5*Val(n_l-2,j+1)+Val(n_l-3,j+1))/48 ...
                                 + ThetPartial;

                p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
                Val(:,j)  = p(1:K);
            end

            %%% Now to Gamm+1
            j = Gamm+1;
            Thet(1:n_l-1)   = 0;
            Thet(n_l)       = ThetPartial;

            p    = ifft(toepM.*fft([Thet; zeros(K,1)]));
            Val(:,j)  = p(1:K);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% CUMULATIVE PARISIAN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else 
        for m=M-2:-1:0

            for j=1:Gamm
                Thet(1)      = (13*Val(1,j+1)+15*Val(2,j+1)-5*Val(3,j+1)+Val(4,j+1))/48;
                Thet(2:n_l-1) = (Val(1:n_l-2,j+1)+10*Val(2:n_l-1,j+1)+Val(3:n_l,j+1))/12;

                Thet(n_l)     = (13*Val(n_l,j+1)+15*Val(n_l-1,j+1)-5*Val(n_l-2,j+1)+Val(n_l-3,j+1))/48 ...
                              + (13*Val(n_l,j)+15*Val(n_l+1,j)-5*Val(n_l+2,j)+Val(n_l+3,j))/48;                          

                Thet(n_l+1:K-1) = (Val(n_l:K-2,j)+10*Val(n_l+1:K-1,j)+Val(n_l+2:K,j))/12;
                Thet(K)        = (13*Val(K,j)+15*Val(K-1,j)-5*Val(K-2,j)+Val(K-3,j))/48;

                p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
                Val(:,j)  = p(1:K);
            end

            %%% Now to Gamm+1
            j = Gamm+1;
            Thet(1:n_l-1)   = 0;
            Thet(n_l)       = (13*Val(n_l,j)+15*Val(n_l+1,j)-5*Val(n_l+2,j)+Val(n_l+3,j))/48;
            Thet(n_l+1:K-1) = (Val(n_l:K-2,j)+10*Val(n_l+1:K-1,j)+Val(n_l+2:K,j))/12;

            p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
            Val(:,j)  = p(1:K);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif down ~=1 && call == 1  %UOC
    %u    = log(H/S_0);
    n_u = n_h;

    sigma = 1 - zeta; sigma_plus = (q_plus-.5)*sigma; sigma_minus = (q_minus-.5)*sigma;

    es1 = exp(dx*sigma_plus); es2 = exp(dx*sigma_minus);

    dbar_0 = .5 + zeta*(.5*zeta-1);
    dbar_1 = sigma*(1 - .5*sigma);

    d_0 = exp((rho+dx)*.5)*sigma^2/18*(5*((1-q_minus)*es2 +(1-q_plus)*es1) + 4);
    d_1 = exp((rho-dx)*.5)*sigma/18*(5*( (.5*(zeta+1) +sigma_minus)*es2 + (.5*(zeta+1) +sigma_plus)*es1 ) + 4*(zeta+1) );
           

    %%%%Thet(1)        =  W/2 - H*varthet_01;
    Thet = zeros(K,1);
    Thet(nbar)     =  W*(exp(-rho)*d_0 - dbar_0);
    Thet(nbar + 1) =  W*(exp(dx-rho)*(varthet_01 +d_1) -(.5 + dbar_1) );
    Thet(nbar + 2:K)  = exp(xmin +dx*(nbar+1:K-1))*S_0*varthet_star - W;

    Val = zeros(K,Gamm+1);
    p   = ifft(toepM.*fft([Thet;zeros(K,1)]));
    for j=1:Gamm
       Val(:,j) = p(1:K); 
    end
    
    Thet(n_u+1:K) = 0;  %after the barrier u, the coefficients are zero
    Thet(n_u) = H*varthet_m10 - .5*W;  %because this is the case where you can knock out at termination
    p   = ifft(toepM.*fft([Thet;zeros(K,1)]));
    Val(:,Gamm+1) = p(1:K);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% RESETTING PARISIAN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% For the resetting, its more efficient to pull the
    %%% Thet(1:n_u-1) ouside of the loop through j = 1:Gamm
    if resetting == 1
        for m=M-2:-1:0
            %NOTE: these must be defined outside of loop, else we change them and then the next function requires the unchanged value!
            
            Thet(1)      = (13*Val(1,1)+15*Val(2,1)-5*Val(3,1)+Val(4,1))/48;
            Thet(2:n_u-1) = (Val(1:n_u-2,1)+10*Val(2:n_u-1,1)+Val(3:n_u,1))/12;
            ThetPartial = (13*Val(n_u,1)+15*Val(n_u-1,1)-5*Val(n_u-2,1)+Val(n_u-3,1))/48 ;
            
            for j=1:Gamm 
                Thet(n_u)     = ThetPartial ...
                    + (13*Val(n_u,j+1)+15*Val(n_u+1,j+1)-5*Val(n_u+2,j+1)+Val(n_u+3,j+1))/48;
                Thet(n_u+1:K-1) = (Val(n_u:K-2,j+1)+10*Val(n_u+1:K-1,j+1)+Val(n_u+2:K,j+1))/12;
                Thet(K)        = (13*Val(K,j+1)+15*Val(K-1,j+1)-5*Val(K-2,j+1)+Val(K-3,j+1))/48;

                p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
                Val(:,j)  = p(1:K);
            end

            %%% Now to Gamm+1
            j = Gamm+1;
            Thet(n_u)     = ThetPartial;
            Thet(n_u+1:K) = 0;

            p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
            Val(:,j)  = p(1:K);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% CUMULATIVE PARISIAN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        for m=M-2:-1:0
            for j=1:Gamm
                Thet(1)      = (13*Val(1,j)+15*Val(2,j)-5*Val(3,j)+Val(4,j))/48;
                Thet(2:n_u-1) = (Val(1:n_u-2,j)+10*Val(2:n_u-1,j)+Val(3:n_u,j))/12;

                Thet(n_u)     = (13*Val(n_u,j)+15*Val(n_u-1,j)-5*Val(n_u-2,j)+Val(n_u-3,j))/48 ...
                              + (13*Val(n_u,j+1)+15*Val(n_u+1,j+1)-5*Val(n_u+2,j+1)+Val(n_u+3,j+1))/48;

                Thet(n_u+1:K-1) = (Val(n_u:K-2,j+1)+10*Val(n_u+1:K-1,j+1)+Val(n_u+2:K,j+1))/12;
                Thet(K)        = (13*Val(K,j+1)+15*Val(K-1,j+1)-5*Val(K-2,j+1)+Val(K-3,j+1))/48;

                p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
                Val(:,j)  = p(1:K);
            end

            %%% Now to Gamm+1
            j = Gamm+1;
            Thet(1)       = (13*Val(1,j)+15*Val(2,j)-5*Val(3,j)+Val(4,j))/48;
            Thet(2:n_u-1) = (Val(1:n_u-2,j)+10*Val(2:n_u-1,j)+Val(3:n_u,j))/12;
            Thet(n_u)     = (13*Val(n_u,j)+15*Val(n_u-1,j)-5*Val(n_u-2,j)+Val(n_u-3,j))/48;
            Thet(n_u+1:K)  =0;

            p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
            Val(:,j)  = p(1:K);
        end
    end
end


if interp_Atend ~= 1
    price = Val(nnot,gamm0);

else  %%% INTERPOLATION
%     Use 5 Point Cubic Interpolation
%     xnot = xmin +(nnot-1)*dx;
%     xs = [xnot-2*dx,xnot-dx,xnot,xnot+dx,xnot+2*dx];
%     ys = [Val(nnot-2,1),Val(nnot-1,1),Val(nnot,1),Val(nnot+1,1),Val(nnot+2,1)];
%     price = spline(xs,ys,0);

    dd = 0 - (xmin+ (nnot -1)*dx); 
    price = Val(nnot,gamm0) + (Val(nnot+1,gamm0) - Val(nnot,gamm0))*dd/dx;
end

end

