function alpha = getTruncationAlpha(T, L1, modelInput, model)
% model: Levy models (BSM, CGMY, NIG, MJD, Kou)
%        Affine models (Heston)

if model == 6  % Heston, already incorpartes T
    alpha = L1*sqrt(abs(modelInput.c2) + sqrt(abs(modelInput.c4)));    % Note: this includes T
else
    alpha = L1*sqrt(abs(modelInput.c2*T) + sqrt(abs(modelInput.c4*T)));    % In this case, we now incorporate T
end
end


