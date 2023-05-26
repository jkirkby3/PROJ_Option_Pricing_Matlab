function [prices, stdErrs] = Price_MC_Return_Barrier_Strikes_func(Spath, call, l, u, Kvec, M, mult, rebate, r, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Calculates Return Barrier option prices for vector of strikes, given the simulatd paths 
%          This version Allows for a rebate. 
%          Note: to price knock-in options, use parity
% Returns: prices and standard errors for each of the supplied strikes
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% Spath = paths of underlying, dimension N_sim x M+1, where M = number of time steps (since includes S_0)
% call = 1 for call (else put)
% l = lower return barrier, e.g. -0.05
% u = upper return barrier, e.g. 0.05
% Kvec = strike vector
% M = number of monitoring points, e.g. 252 for "daily" monitoring
% mult = time partitioning multiplier to reduce bias (e.g. mult = 1, 2,..)
% rebate = rebate which is paid upon barrier breach, e.g. 5.0
% r = Interest rate
% q = dividend yield
% T = Time (in years)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prices = zeros(length(Kvec),1);
stdErrs = zeros(length(Kvec),1);
N_sim = size(Spath,1);  % number of paths

knock_time = zeros(N_sim,1);  % time of knock-out

M_mult = M*mult;  %time partitioning to reduce bias
dt_mult = T/M_mult;

for n = 1:N_sim
    Slast = Spath(n,1);
    for m=1+mult:mult:M_mult+1   % We are calculating returns, so start at the next monitoring time
        Snext = Spath(n,m);
        ret = log(Snext/Slast);
        if ret < l || ret > u
            knock_time(n) = (m-1)*dt_mult;
            break
        end
        Slast = Snext;
    end
end

if rebate > 0
   disc_rebate = rebate*exp(-r*knock_time).*(knock_time>0);   % discounted rebate
else
   disc_rebate = 0; 
end

for k = 1:length(Kvec)
    K = Kvec(k);
    if call ==1
    	payoffs  = exp(-r*T)*max(0, Spath(:,M_mult+1) - K).*(knock_time==0) + disc_rebate;
        
    else
        payoffs  = exp(-r*T)*max(0, K - Spath(:,M_mult+1)).*(knock_time==0) + disc_rebate;
    end
    
    prices(k) = mean(payoffs);
    stdErrs(k) = std(payoffs) / sqrt(N_sim);
end


end

