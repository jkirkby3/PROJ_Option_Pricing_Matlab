function chf = rnCHF_disc(z, levyExponentRN, H, G, T, n, stateGrid, RNdrift)
    % k_0 is the bracketing index: grid(k_0) <= v0 < grid(k_0+1)
    %Returns the chf with is linear interpolation of the two bracketing chfs
    dt = T/n;   % different than T...
    P = expm(dt*G);
    
    numStates = length(stateGrid);
    numZ = length(z);
    
    ones_ = ones(numStates, 1);   % column vector
    chf = zeros(numStates, numZ); % chf for each initial state 

    for j = 1:numZ   % populate one column (one theta(i) value) at a time
        E = diag(exp(levyExponentRN(z(j))*H*dt));
        chf(:, j) = E*ones_ * exp(1i*RNdrift*z(j));   % to initialize matrix vector mult
        EP = E*P;
        for k = 1:n
            chf(:, j) = EP*chf(:, j); 
        end
     
    end
    
end