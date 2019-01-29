function chf = rnCHF_cont(z, levyExponentRN, H, G, T, stateGrid, RNdrift)
    numStates = length(stateGrid);
    numZ = length(z);
    
    ones_ = ones(numStates, 1);   % column vector
    chf = zeros(numStates, numZ); % chf for each initial state 
    
    for j = 1:numZ   % populate one column (one theta(i) value) at a time
        chf(:, j) = expm( T*(G + diag(levyExponentRN(z(j))*H))) * ones_ * exp(1i*RNdrift*z(j));   % to initialize matrix vector mult
    end
    
end