function price = DiscreteVariance_PROJ_Levy( N,alph,M,r,T,K,rnCHF,contract )
% N = #basis points
% alph = log-asset grid width param
% M = # Monitoring dates (not including S_0)
% r = interest rate
% T = time to maturity
% K = strike (only matters for an option, but is always required)

dx   = 2*alph/(N-1);
a    = 1/dx;
A    = 32*a^4;
C_aN = A/N;
xmin = (1-N/2)*dx; 

dxi    = 2*pi*a/N;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% PSI Matrix: 5-Point GAUSSIAN
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

NNM = N;  %% in this case, same dimension

PSI     = zeros(N,NNM);    %The first row will remain ones (since the chf will always be one at that point)
PSI(1,:) = ones(1,NNM);

%%%% Sample
Neta  = 5*(NNM) + 15;   %sample size
Neta5 = (NNM) + 3;
g2    = sqrt(5-2*sqrt(10/7))/6;
g3    = sqrt(5+2*sqrt(10/7))/6;
v1    = .5*128/225; v2 = .5*(322+13*sqrt(70))/900;  v3 = .5*(322 - 13*sqrt(70))/900;

thet                 = zeros(1,Neta);   %sample initialized
thet(5*(1:Neta5)-2)  = xmin -1.5*dx + dx*(0:Neta5-1);
thet(5*(1:Neta5)-4)  = xmin -1.5*dx + dx*(0:Neta5-1) - dx*g3;
thet(5*(1:Neta5)-3)  = xmin -1.5*dx + dx*(0:Neta5-1) - dx*g2;
thet(5*(1:Neta5)-1)  = xmin -1.5*dx + dx*(0:Neta5-1) + dx*g2;
thet(5*(1:Neta5))    = xmin -1.5*dx + dx*(0:Neta5-1) + dx*g3;

%%%% Weights
sig      = [-1.5-g3, -1.5-g2, -1.5, -1.5+g2, -1.5+g3, -.5-g3, -.5-g2, -.5, -.5+g2, -.5+g3,];
sig(1:5) = (sig(1:5) + 2).^3/6;
sig(6:10) = 2/3 - .5*(sig(6:10)).^3 - (sig(6:10)).^2;

sig([1 5 6 10]) = v3*sig([1 5 6 10]); sig([2 4 7 9]) = v2*sig([2 4 7 9]); sig([3 8]) = v1*sig([3 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%NEW STEP: multiple sig by Upsilon_{a,N}
sig = C_aN*sig;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fill Matrix
%%%% NOTE: this can be made MORE EFFICIENT by using symmetery of x^2

zz  = exp(1i*dxi*thet.^2); %% in general, 1i*dxh(thet)
thet   = zz; 

for j=2:N
    PSI(j,:) =  sig(1)*(thet(1:5:Neta-19) + thet(20:5:Neta)) ...
              + sig(2)*(thet(2:5:Neta-18) + thet(19:5:Neta-1)) ...
              + sig(3)*(thet(3:5:Neta-17)  + thet(18:5:Neta-2)) ...
              + sig(4)*(thet(4:5:Neta-16)  + thet( 17:5:Neta-3)) ...
              + sig(5)*(thet(5:5:Neta-15)  + thet( 16:5:Neta-4)) ...
              + sig(6)*(thet(6:5:Neta-14)  + thet( 15:5:Neta-5)) ...
              + sig(7)*(thet(7:5:Neta-13)  + thet( 14:5:Neta-6)) ...
              + sig(8)*(thet(8:5:Neta-12)  + thet( 13:5:Neta-7)) ...
              + sig(9)*(thet(9:5:Neta-11)  + thet( 12:5:Neta-8)) ...
              + sig(10)*(thet(10:5:Neta-10)  + thet( 11:5:Neta-9));

    thet = thet.*zz;
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Find phi_{Y_1}
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

xi    = dxi*(1:N-1)';  %REDEFINED from above, excludes the zero at xi_1

b0    = 1208/2520; b1 = 1191/2520; b2 = 120/2520; b3 = 1/2520;
zeta  = (sin(xi/(2*a))./xi).^4./(b0 + b1*cos(xi/a) +b2*cos(2*xi/a) +b3*cos(3*xi/a));
hvec  = exp(-1i*xmin*xi).*zeta;  %dont define first element


AA = 1/A;
beta  = [AA; rnCHF(xi).*hvec];   %grand(end)=.5*grand(end);
beta  = real(fft(beta));

phi = PSI(:,1:N)*beta;  

%%%% FIND
phi = phi.^M;  %NOTE: this is the entire convolution step, b/c indepdendent stationary increments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Redfine xmin for the final inversion

if contract == 1 || contract == 3 %Variance Swap or Call on Variance   
    grid = dx*(0:N-1);  %NOTE: this holds for call on variance too, since we subtract out the KT
    grid(1) = dx/6;
end

if contract == 1 %Variance Swap
    xmin = 0;
elseif contract == 3 %Call on Variance
    xmin = K*T;
end

C_aN  = 24*a^2/N;
dw    = 2*pi*a/N;
grand = dw*(1:N-1)';
grand = exp(-1i*xmin*grand).*phi(2:N).*(sin(grand/(2*a))./grand).^2./(2+cos(grand/a));

%%% Test with FILTER
applyFilter = 0;
if applyFilter == 1
    epsM = 1.2204e-16;   %matlabs machine epsilon
    alphaeps = -log(epsM);
    pp = 8; %order of the filter, must be even
    filter = exp(-alphaeps*(xi/(2*pi*a)).^pp);
    grand = grand.*filter;
end

beta  = real(fft([1/(24*a^2); grand]));   %%%%  NOTE: all toep matrices incorporate exp(-r*dt)


if contract == 1 || contract == 3 %Variance Swap or call
    price = grid(1:N/2)*beta(1:N/2);
    price = C_aN*price;
end
    
if contract == 1 %Variance Swap
    price = price/T;  %can generalize to include anualization factor if needed
elseif contract == 3 %Variance Call
    price = exp(-r*T)*price/T;
else
    fprintf('Only contract types 1 and 3 are currently supported \n')
end

end

