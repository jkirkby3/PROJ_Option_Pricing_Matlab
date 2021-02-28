function Spath = Simulate_RegimeSwitching_Diffusion_func( N_sim, M, T, S_0, drift_vec, sigma_vec, Q, initial_state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Simulates Paths of Regime Switching Diffusion Models
%        Uses log-Euler scheme
% Returns: paths of dimension (N_sim, M+1), since they include S_0 
%          ... Simulates N_sim paths, each row is a full path starting from S_0, ending with S_M (M+1 points in path)
% Author: Justin Lars Kirkby
%
% -----------------
% Params
% -----------------
% N_sim = # paths
% M = #time steps on [0,T], time step is dt=T/M, so each path has M+1 points
% T = time to maturity, ie path is on [0,T]
% S_0 = initial underlying value (e.g. S_0=100)
% drift_vec = vector of drift coefficient by regime state, e.g. r_i - q_i, where r is interest rate in state i, and q is div yield
% sigma_vec = diffusion coefficients in each state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = T/M;

% Initialize the CDFs used to simulate transitions between regimes
cdfs = cumsum(expm(Q*dt).', 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

Spath = zeros(N_sim,M+1);
Spath(:,1) = S_0;
Sigsqdt_vec = reshape(sigma_vec*sqrt(dt), length(Q),1);
drift_vec = reshape((drift_vec - .5*sigma_vec.^2)*dt,  length(Q),1);  % Compensate the drift

prev_states = initial_state*ones(N_sim,1);


for m = 1:M
    W1 = randn(N_sim,1); 

    % Simulate the next Regime
    states = sim_ctmc_step(prev_states, cdfs);
    prev_states = states;

    Sigsqdt = Sigsqdt_vec(states);
    drift = drift_vec(states);

    % Simulate Underlying With this new regime
    Spath(:,m+1) = Spath(:,m).*exp(drift + Sigsqdt.*W1);  %log scheme

end
end

function k = sim_ctmc_step(states, cdfs)
    u = rand(size(states));
    [~,k] =max(u <= cdfs(states,:),[], 2);   % max here will find the FIRST among those which satisfy the condition
end


%%%%%%%%%%%%%%%%%%
% Non-Vectorized Version (Slow)
%%%%%%%%%%%%%%%%%%

% function k = sim_ctmc_step(i, cdfs)
%     % TODO: use find instead   j_0 = find(v >= v0, 1, 'first');
%     u = rand;
%     k =find(u <= cdfs(i,:), 1, 'first');
% end

% for m = 1:M
%     W1 = randn(N_sim,1); 
%     for j = 1: N_sim
%         % Simulate the next Regime
%         state = sim_ctmc_step(prev_states(j), cdfs);
%         prev_states(j) = state;
%         
%         Sigsqdt = Sigsqdt_vec(state);
%         drift = drift_vec(state);
%         
%         % Simulate Underlying With this new regime
%         Spath(j,m+1) = Spath(j,m).*exp(drift + Sigsqdt*W1(j));  %log scheme
%     end
% end







