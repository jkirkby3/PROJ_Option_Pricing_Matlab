%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: This script and pricing function are in progress, and have not been well tested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating forward value of the call with Antonov's mapping strategy
T      = 1;  %Time to maturity
r      = 0.00;    %Risk-free interest rate
F_0    = 1.1;    %Initial forward value

ModParams.v0     = 0.2;  %Inital volatility
ModParams.beta   = 0.7;  %Exponent
ModParams.alpha  = 0.08;  %Vol-vol
ModParams.rho    = 0; %Correlation

% ModParams.v0     = 0.25;  %Inital volatility
% ModParams.beta   = 0.6;  %Exponent
% ModParams.alpha     = 0.3;  %Vol-vol
% ModParams.rho     = -0.5; %Correlation


Kvec   = F_0*[0.6 0.8 0.90 0.95 0.999 1.05 1.10 1.2 1.4];
call   = 0;

%%% Price Strikes
prices = zeros(length(Kvec),1);
for k=1:length(Kvec)
    K = Kvec(k);
    prices(k) = SABR_European_AntonovApprox(F_0,K,T,call,r,ModParams);
    fprintf('%.8f\n', prices(k));
end
