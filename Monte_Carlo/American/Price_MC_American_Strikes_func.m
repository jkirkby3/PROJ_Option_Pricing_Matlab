function [prices, stdErrs] = Price_MC_American_Strikes_func(Spath, disc, call, Kvec, polyOrder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Longstaff-Schwartz Algo: Calculates American option prices for vector of strikes, given the simulatd paths 
% Returns: prices and standard errors for each of the supplied strikes
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% Spath = paths of underlying, dimension N_sim x M+1, where M = number of time steps (since includes S_0)
% disc = discount factor for each time step, e.g. exp(-r*dt), where dt is time step, and r is interest rate
% call = 1 for call (else put)
% Kvec = strike vector
% polyOrder = order of polynomial regression of continuation val, e.g. 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    polyOrder = 3;  % Cubic polynomial by default
end

prices = zeros(length(Kvec), 1);
stdErrs = zeros(length(Kvec), 1);

for k = 1:length(Kvec)
    K = Kvec(k);
    [prices(k), stdErrs(k)]  = Price_MC_American_func(Spath, disc, call, K, polyOrder);
end

end

function [price, stdErr] = Price_MC_American_func(Spath, disc, call, K, polyOrder)
    M = size(Spath, 2) - 1;  % number of time steps 
    N_sim = size(Spath, 1);   % number of paths
    
    if call == 1
    	payoff = max(Spath(:,M+1) - K,0);
    else
        payoff = max(K-Spath(:,M+1),0);
    end
    
    for m=M:-1:2
        payoff = payoff*disc; %holding value
        
        if call == 1
            EV = max(Spath(:,m) - K, 0); %exercise value at tn
        else
            EV = max(K - Spath(:,m), 0); %exercise value at tn
        end
        
        index = find(EV>0);  % just regress on in the money positions, otw there is cluster at zero due to payoff

        regression = polyfit(Spath(index, m), payoff(index), polyOrder);
        EHV = polyval(regression, Spath(index, m));

        si = size(index);

        for j=1:si(1)
            index_j = index(j);
            if EHV(j) < EV(index_j) 
                payoff(index_j) = EV(index_j);
            end
        end
    end
        
    price = mean(payoff)*disc;
    stdErr = disc*std(payoff) / sqrt(N_sim);
end
