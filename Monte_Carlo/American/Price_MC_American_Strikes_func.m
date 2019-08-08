function [prices, stdErrs] = Price_MC_American_Strikes_func(Spath, disc, call, Kvec )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Longstaff-Schwartz Algo: Calculates American option prices for vector of strikes, given the simulatd paths 
% Returns: prices and standard errors for each of the supplied strikes
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% disc = discount factor, e.g. exp(-r*T)
% call = 1 for call (else put)
% Kvec = strike vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prices = zeros(length(Kvec), 1);
stdErrs = zeros(length(Kvec), 1);

for k = 1:length(Kvec)
    K = Kvec(k);
    [prices(k), stdErrs(k)]  = Price_MC_American_func(Spath, disc, call, K);
end

end

function [price, stdErr] = Price_MC_American_func(Spath, disc, call, K )
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
        
        index = find(EV>0);

        regression = polyfit(Spath(index,m), payoff(index), 4);
        EHV = polyval(regression, Spath(index,m));

        si = size(index);

        for j=1:si(1)
            if EHV(j) < EV(index(j)) 
                payoff(index(j)) = EV(index(j));
            end
        end
    end
        
    price = mean(payoff)*disc;
    stdErr = disc*std(payoff) / sqrt(N_sim);
end
