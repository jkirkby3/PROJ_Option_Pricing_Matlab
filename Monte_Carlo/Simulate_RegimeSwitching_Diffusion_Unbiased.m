function Svals = Simulate_RegimeSwitching_Diffusion_Unbiased( N_sim, T, S_0, drift_vec, sigma_vec, Q, initial_state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Simulates Terminal Value (S_T) of Regime Switching Diffusion Models
%        This scheme is unbiased, but only generates the terminal S values
% Returns: vector of size N_sim
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% N_sim = # samples
% T = time to maturity, ie path is on [0,T]
% S_0 = initial underlying value (e.g. S_0=100)
% drift_vec = vector of drift coefficient by regime state, e.g. r_i - q_i, where r is interest rate in state i, and q is div yield
% sigma_vec = diffusion coefficients in each state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Svals = zeros(N_sim,1);
lambdas = -diag(Q);  % Transition rates
P = (Q + diag(lambdas))./lambdas;  % Remove diagonal, no self-transition 
cdfs = cumsum(P, 2);

for n = 1:N_sim
   S = S_0;
   t_last = 0;
   state = initial_state;
   while t_last < T
       tau = exprnd(lambdas(state));
       t = min(T, t_last + tau);
       tau = t - t_last;
       
       Sigsqdt = sigma_vec(state)*sqrt(tau);
       drift = (drift_vec(state) - 0.5*sigma_vec(state)^2)*tau;
       
       S = S.*exp(drift + Sigsqdt*randn());  %log scheme
       state = next_state(cdfs(state, :));
       t_last = t;
   end
   Svals(n) = S;
end

end

function j = next_state(cdf)
u = rand();
j = find(u <= cdf, 1);
end
