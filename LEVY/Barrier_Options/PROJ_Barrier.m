function price = PROJ_Barrier(N, alph, call, down, S_0, W, H, M, r, q, rnCHF, T, rebate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Pricing Function for Discrete Barrier Options using PROJ method
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
% Author: Justin Lars Kirkby
%
% ----------------------
% Contract/Model Params 
% ----------------------
% S_0 = initial stock price (e.g. 100)
% W   = strike  (e.g. 100)
% r   = interest rate (e.g. 0.05)
% q   = dividend yield (e.g. 0.05)
% T   = time remaining until maturity (in years, e.g. T=1)
% M   = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
% call = 1 for call (else put)
% down = 1 for down and out (otherwise it's up and out)
% H    = barrier
% rebate = rebate paid immediately upon passing the barrier (knocking-out)
% rnCHF = risk netural characteristic function (function handle with single argument)
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% alph  = grid with is 2*alph
% N     = number of grid/basis points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%*****************************
%TODO: Refactor similar to parisian code to simplify
%*****************************

if nargin < 13
    rebate = 0;
end

mult = 1;  % How much overhang to remove aliasing (hardcoded to one for single barriers)
interp_Atend = 0;  %see grid determination: this may get set to 1 below, indicating an interpolation at end between V(nnot) and V(nnot+1)


K = N/2;

dx = 2*alph/(N-1); a = 1/dx;

dt   = T/M;
nrdt = -r*dt;  nqdt = -q*dt;

Thet  = zeros(K,1);

%%%%%%  Gaussian Quad Constants
q_plus = (1 + sqrt(3/5))/2;  q_minus = (1 - sqrt(3/5))/2;
b3  = sqrt(15); b4 = b3/10;

if down ==1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%  DOWN & OUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    l    = log(H/S_0);
    xmin = l;
    nnot =  floor(1-xmin*a);
    
    if nnot>=K 
        fprintf('nnot is %.0f while K is %.0f, need to increase alpha \n',nnot, K)
    end
    
    if call == 1 && nnot == 1 %In this case a DOC with H near S_0 is still very valuable, so setting alph too small is bad idea        
        interp_Atend = 1; %Instruct to use interpolation at algorithm end
        %no change is made to dx
    else
        nnot = max(2, floor(1-xmin*a)); dx = l/(1-nnot);
    end    
    a    = 1/dx;

    lws   = log(W/S_0);
    nbar  = floor(a*(lws - xmin)+1);
    rho   = lws - (xmin+(nbar - 1)*dx);
    zeta  = a*rho;
    
    a2    = a^2;
    zmin  = (1 - K)*dx;  %Kbar corresponds to zero
   
    %%% Extend Pbar, only to invert (aliasing)
    Nmult = mult*N;
    Cons = 24*a2*exp(nrdt)/Nmult;
    dw    = 2*pi*a/Nmult;
    grand = dw*(1:Nmult-1);
    grand = exp(-1i*zmin*grand).*rnCHF(grand).*(sin(grand/(2*a))./grand).^2./(2+cos(grand/a));
    beta  = Cons*real(fft([1/(24*a2) grand]));  
    
    toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];   toepM = fft(toepM);
    
    %%%% PAYOFF CONSTANTS-----------------------------------
    varthet_01 = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
    varthet_m10 = exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4)/18;
    varthet_star = varthet_01 + varthet_m10;
    %varthet_star = (cosh(dx/2)*(4 + 5*cosh(b4*dx) - b3*sinh(dx/2)*sinh(b4*dx)))/9;
    %%%%---------------------------------------------------
    
    if rebate ~= 0
        val_rebate = rebate*[ fliplr(cumsum(beta(1:1:K-1)))';0];  % NOTE: this includes the discounting via beta
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %%%% DOC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if call == 1  
        sigma = 1 - zeta; sigma_plus = (q_plus-.5)*sigma; sigma_minus = (q_minus-.5)*sigma;

        es1 = exp(dx*sigma_plus); es2 = exp(dx*sigma_minus);

        dbar_0 = .5 + zeta*(.5*zeta-1);
        dbar_1 = sigma*(1 - .5*sigma);

        d_0 = exp((rho+dx)*.5)*sigma^2/18*(5*((1-q_minus)*es2 +(1-q_plus)*es1) + 4);
        d_1 = exp((rho-dx)*.5)*sigma/18*(5*( (.5*(zeta+1) +sigma_minus)*es2 + (.5*(zeta+1) +sigma_plus)*es1 ) + 4*(zeta+1) );

        
        Thet(nbar)          = W*(exp(-rho)*d_0 - dbar_0);
        Thet(nbar + 1)      = W*(exp(dx-rho)*(varthet_01 +d_1) -(.5 + dbar_1) );
        Thet(nbar + 2:K)    = exp(xmin +dx*(nbar+1:K-1))*S_0*varthet_star - W;
        
        Thet(1) = Thet(1) + 0.5*rebate;
        
        toepR = [beta(2*K:-1:K +1)'; 0 ;zeros(K-1,1)];   toepR = fft(toepR);
        
        %%% In the next two lines, exp(r*dt) is to compensate for
        %%% multiplication by Cons2 in each of the toep matrices
        Thetbar1 = exp(r*dt)*W*cumsum(beta(2*K:-1:K +1))';
        Thetbar2 = exp(r*dt)*S_0*varthet_star*exp(xmin + dx*(K:2*K -1))';
        p = ifft(toepR.*fft([Thetbar2; zeros(K,1)]));
        Thetbar2 = p(1:K);

        %%%%%%%
        p = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));

        if rebate ~= 0  
             Val = p(1:K) + exp(nrdt)*(Thetbar2-Thetbar1)+ val_rebate;
        else
             Val = p(1:K) + exp(nrdt)*(Thetbar2-Thetbar1);
        end
        %%%%%%%

        for m=M-2:-1:0
            Thet(2:K -1) = (Val(1:K-2)+10*Val(2:K-1)+Val(3:K))/12;
            Thet(1)      = (13*Val(1)+15*Val(2)-5*Val(3)+Val(4))/48;
            Thet(K)      = 2*(13*Val(K)+15*Val(K-1)-5*Val(K-2)+Val(K-3))/48;  %NOTE: 2*theta(K) b/c of augmenting
            
            Thet(1) = Thet(1) + 0.5*rebate;   % account for overhang into the knock-out region
            %%%%%%%
            p        = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
            Val(1:K) = p(1:K)+ exp(nqdt*(M - m -1))*Thetbar2 - exp(nrdt*(M - m -1))*Thetbar1;
            
            if rebate ~= 0
                Val = Val + val_rebate;
            end
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %%%% DOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else 
        zeta_plus = zeta*q_plus; zeta_minus = zeta*q_minus;
        rho_plus = rho*q_plus; rho_minus = rho*q_minus;
        
        ed1 = exp(rho_minus); ed2 = exp(rho/2); ed3 = exp(rho_plus);

        dbar_1 = zeta^2/2;
        dbar_0 = zeta - dbar_1;         %  dbar_1 = zeta + .5*((zeta - 1)^2 - 1);
        d_0    = zeta*(5*( (1-zeta_minus)*ed1 + (1-zeta_plus)*ed3 ) + 4*(2-zeta)*ed2)/18;
        d_1    = zeta*( 5*(zeta_minus*ed1 + zeta_plus*ed3) + 4*zeta*ed2 )/18;            

        Thet(1)        =  W/2 - H*varthet_01;
        Thet(2:nbar-1) =  W - exp(xmin +dx*(1:nbar-2))*S_0*varthet_star;
        Thet(nbar)     =  W*(.5 + dbar_0 - exp(-rho)*(varthet_m10 + d_0));
        Thet(nbar + 1) =  W*(dbar_1 - exp(- rho)*d_1);
        
        Thet(1) = Thet(1) + 0.5*rebate;
        %%%%%%%
        p   = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
        if rebate ~= 0  
             Val = p(1:K) + val_rebate;
        else
             Val = p(1:K);
        end
        %%%%%%%

        for m=M-2:-1:0
            Thet(1)      = (13*Val(1)+15*Val(2)-5*Val(3)+Val(4))/48;
            Thet(K)      = (13*Val(K)+15*Val(K-1)-5*Val(K-2)+Val(K-3))/48;
            Thet(2:K -1) = (Val(1:K-2)+10*Val(2:K-1)+Val(3:K))/12;
            
            Thet(1) = Thet(1) + 0.5*rebate;   % account for overhang into the knock-out region
            %%%%%%%
            p    = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));
            Val  = p(1:K);
            
            if rebate ~= 0
                Val = Val + val_rebate;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%  UP & OUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else 
    lws  = log(W/S_0);
    u    = log(H/S_0);
    %nnot = min(K-1, floor(K-a*u)); %The min will always hold
    nnot = floor(K-a*u);
    if call ~= 1 && nnot == K-1
        interp_Atend= 1; %dont change value of dx
    else
        dx   = u/(K-nnot);  a    = 1/dx;
    end   
    
    xmin = u - (K-1)*dx;  %%%%NOTE: this used to be  xmin = u - (K-2)*dx;
    nbar = floor(a*(lws - xmin)+1);
    
    rho   = lws - (xmin+(nbar - 1)*dx);
    zeta  = a*rho;
    
    a2    = a^2;
    zmin  = (1 - K)*dx;  %Kbar corresponds to zero
   
    %%% Extend Pbar, only to invert (aliasing)
    Nmult = mult*N;
    Cons = 24*a2*exp(nrdt)/Nmult;
    dw    = 2*pi*a/Nmult;
    grand = dw*(1:Nmult-1);
    grand = exp(-1i*zmin*grand).*rnCHF(grand).*(sin(grand/(2*a))./grand).^2./(2+cos(grand/a));
    beta  = Cons*real(fft([1/(24*a2) grand]));  
    
    toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];
    toepM = fft(toepM);

    %%%% PAYOFF CONSTANTS-----------------------------------
    varthet_01 = exp(.5*dx)*(5*cosh(b4*dx) - b3*sinh(b4*dx) + 4)/18;
    varthet_m10 = exp(-.5*dx)*(5*cosh(b4*dx) + b3*sinh(b4*dx) + 4)/18;
    varthet_star = varthet_01 + varthet_m10;
    %varthet_star = (cosh(dx/2)*(4 + 5*cosh(b4*dx) - b3*sinh(dx/2)*sinh(b4*dx)))/9;
    
    if rebate ~= 0
        val_rebate = rebate*cumsum(beta(2*K:-1:K +1))';  % NOTE: this includes the discounting via beta
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %%%% UOC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if call == 1 
        
        sigma = 1 - zeta; sigma_plus = (q_plus-.5)*sigma; sigma_minus = (q_minus-.5)*sigma;
        es1 = exp(dx*sigma_plus); es2 = exp(dx*sigma_minus);

        dbar_0 = .5 + zeta*(.5*zeta-1);
        dbar_1 = sigma*(1 - .5*sigma);

        d_0 = exp((rho+dx)*.5)*sigma^2/18*(5*((1-q_minus)*es2 +(1-q_plus)*es1) + 4);
        d_1 = exp((rho-dx)*.5)*sigma/18*(5*( (.5*(zeta+1) +sigma_minus)*es2 + (.5*(zeta+1) +sigma_plus)*es1 ) + 4*(zeta+1) );

        
        Thet(nbar)          = W*(exp(-rho)*d_0 - dbar_0);
        Thet(nbar + 1)      = W*(exp(dx-rho)*(varthet_01 +d_1) -(.5 + dbar_1) );
        Thet(nbar + 2:K-1)  = exp(xmin +dx*(nbar+1:K-2))*S_0*varthet_star - W;
        Thet(K)             = H*varthet_m10 - .5*W;
        
        Thet(K) = Thet(K) + 0.5*rebate;  % account for overlap of right most basis point with knock-out region
        %%%%%%%
        
        p = ifft(toepM.*fft([Thet(1:K); zeros(K,1)]));
        if rebate ~= 0  
             Val = p(1:K) + val_rebate;
        else
             Val = p(1:K);
        end
        
        %%%%%%%

        for m=M-2:-1:0
            Thet(1)      = (13*Val(1)+15*Val(2)-5*Val(3)+Val(4))/48;
            Thet(K)      = (13*Val(K)+15*Val(K-1)-5*Val(K-2)+Val(K-3))/48;
            Thet(2:K -1) = (Val(1:K-2)+10*Val(2:K-1)+Val(3:K))/12;
            
            Thet(K) = Thet(K) + 0.5*rebate;
            %%%%%%%
            p        = ifft(toepM.*fft([Thet(1:K); zeros(K,1)]));
            Val = p(1:K);
            if rebate ~= 0 
                Val = Val + val_rebate;
            end
        end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %%%% UOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else %(UOP)
        zeta_plus = zeta*q_plus; zeta_minus = zeta*q_minus;
        rho_plus = rho*q_plus; rho_minus = rho*q_minus;
        ed1 = exp(rho_minus); ed2 = exp(rho/2); ed3 = exp(rho_plus);

        dbar_1 = zeta^2/2;
        dbar_0 = zeta - dbar_1;         %  dbar_1 = zeta + .5*((zeta - 1)^2 - 1);
        d_0    = zeta*(5*( (1-zeta_minus)*ed1 + (1-zeta_plus)*ed3 ) + 4*(2-zeta)*ed2)/18;
        d_1    = zeta*( 5*(zeta_minus*ed1 + zeta_plus*ed3) + 4*zeta*ed2 )/18;            

        Thet(1:nbar-1) =  W - exp(xmin +dx*(0:nbar-2))*S_0*varthet_star;
        Thet(nbar)     =  W*(.5 + dbar_0 - exp(-rho)*(varthet_m10 + d_0));
        Thet(nbar + 1) =  W*(dbar_1 - exp(- rho)*d_1);

        Thet(K) = 0.5*rebate;  % account for overlap of right most basis point with knock-out region
        
        toepL = [zeros(K,1); 0 ; beta(K-1:-1:1)'];
        toepL = fft(toepL);

        Thetbar1 = exp(r*dt)*W*[ fliplr(cumsum(beta(1:1:K-1)))';0];
        Thetbar2 = exp(r*dt)*S_0*varthet_star*exp(xmin - dx*(K:-1:1))';
        p = ifft(toepL.*fft([Thetbar2; zeros(K,1)]));
        Thetbar2 = p(1:K);

        p = ifft(toepM.*fft([Thet(1:K); zeros(K,1)]));
        %%%%%%%      
        if rebate ~= 0
            Val = p(1:K) + exp(-r*dt)*(Thetbar1-Thetbar2) + val_rebate;
        else
            Val = p(1:K) + exp(-r*dt)*(Thetbar1-Thetbar2);
        end
        
        %%%%%%%

        %%%%%%%     
        for m=M-2:-1:0
            Thet(1)      = 2*(13*Val(1)+15*Val(2)-5*Val(3)+Val(4))/48;  %NOTE: 2*theta(1) b/c of augmenting
            Thet(K)      = (13*Val(K)+15*Val(K-1)-5*Val(K-2)+Val(K-3))/48;  
            Thet(2:K -1) = (Val(1:K-2)+10*Val(2:K-1)+Val(3:K))/12;
            
            Thet(K) = Thet(K) + 0.5*rebate;   % account for overhang into the knock-out region
            %%%%%%%
            p    = ifft(toepM.*fft([Thet(1:K); zeros(K,1)]));
            Val  = p(1:K) + exp(nrdt*(M - m -1))*Thetbar1 - exp(nqdt*(M - m -1))*Thetbar2 ;
            
            if rebate ~= 0
                Val = Val + val_rebate;
            end
        end

    end
end

if interp_Atend == 1
    dd = 0 - (xmin+ (nnot -1)*dx); 
    price = Val(nnot) + (Val(nnot+1) - Val(nnot))*dd/dx; %ie linear interp of nnot and nnot+1
else
    price = Val(nnot); 
end

end
