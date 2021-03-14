function price = Hilbert_European_Price(h, N, r, q, T, S_0, W, call, rnCHF )
% h = step size
% N = budget
gridL = -N/2:-1; gridR = -fliplr(gridL);
g = @(z) exp(-1i*z*log(W/S_0)).*(S_0.*rnCHF(z-1i) - W*rnCHF(z));
H = sum(g(h*gridL).*(cos(pi*gridL)-1)./gridL + g(h*gridR).*(cos(pi*gridR)-1)./gridR)/pi;

% Call option price
price = .5*real(S_0*exp(-q*T) - W*exp(-r*T)...
        + 1i*exp(-r*T)*H);

if call ~= 1  % price put using put-call parity
    price = price - (S_0*exp(-q*T) - W*exp(-r*T));
end


end

