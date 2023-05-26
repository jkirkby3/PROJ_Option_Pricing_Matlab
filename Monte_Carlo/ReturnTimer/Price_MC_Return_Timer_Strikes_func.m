function [prices, stdErrs] = Price_MC_Return_Timer_Strikes_func(Spath, call, l, u, Kvec, M, mult, r, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Calculates Return Timer option prices for vector of strikes, given the simulatd paths 
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
% r = Interest rate
% q = dividend yield
% T = Time (in years)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prices = zeros(length(Kvec),1);
stdErrs = zeros(length(Kvec),1);
N_sim = size(Spath,1);  % number of paths

knock_time = T*ones(N_sim,1);  % time of knock-out, default to final period
Svals = zeros(N_sim,1);

M_mult = M*mult;  %time partitioning to reduce bias
dt_mult = T/M_mult;

for n = 1:N_sim
    Slast = Spath(n,1);
    knock = 0;
    for m=1+mult:mult:M_mult+1   % We are calculating returns, so start at the next monitoring time
        Snext = Spath(n,m);
        ret = log(Snext/Slast);
        if ret < l || ret > u
            knock = 1;
            knock_time(n) = (m-1)*dt_mult;
            % The return is capped at the time Timer Hits
            if ret < l
                Svals(n) = Slast*exp(l);
            else
                Svals(n) = Slast*exp(u);
            end
            break
        end
        Slast = Snext;
    end
    if knock == 0 % If you dont knock out
        Svals(n) = Spath(n,M_mult+1);
    end
    
end

for k = 1:length(Kvec)
    K = Kvec(k);
    if call ==1
    	payoffs  = exp(-r*knock_time).*max(0, Svals - K);
    else
        payoffs  = exp(-r*knock_time).*max(0, K - Svals);
    end
    
    prices(k) = mean(payoffs);
    stdErrs(k) = std(payoffs) / sqrt(N_sim);
end


end

