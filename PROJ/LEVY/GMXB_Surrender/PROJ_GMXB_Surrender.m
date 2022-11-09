function price = PROJ_GMXB_Surrender(T, M, gmdb_params, modelInputs, N, alph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Price Gauranteed Minimum death benefits with period fees and early surrender in Levy Models
%              using the PROJ method
%
% Author:      Justin Lars Kirkby
% References:  (1) Valuation and optimal surrender of variable annuities
%                   with guaranteed minimum benefits and periodic fees, 
%                   Kirkby and Aguilar 2022, Scandinavian Actuarial Journal
%
%              (2) Efficient Option Pricing By Frame Duality with The Fast
%              Fourier Transform, SIAM J. Financial Math., 2015
%
% Models Supported: Levy Processes, including jump diffusions and Black-Scholes model
% Returns: price of contract
%
% ----------------------
% Contract/Model Params
% ----------------------
% T   = time remaining until maturity
% M   = number of subintervals of [0,T] (total of M+1 monitoring points in time grid, including S_0)
%
% modelInputs
% ------------------------
% r     = interest rate
% q     = dividend yield
% rnCHF = risk netural characteristic function (function handle with single argument)
%
% gmdb_params
% ------------------------
% F_0       = initial fund value
% alpha_fee = period fee rate
% gamma     = surrender penalty rate, 1.0 = 100% fund lost upon surrender, 0.0 = 0% lost, ie no penalty
% g         = floor on growth
% c         = cap on growth
% death_prob_cond = conditional probability of death (mortality) table, based on age
%                   To construct see function make_mortality_table_pmf, set conditional=1
%
% ----------------------
% Numerical (PROJ) Params 
% ----------------------
% alph =  grid with is 2*alph
% N  = number of grid points (power of 2, e.g. 2^12), resolution = 2*alph/(N-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt   = T/M;

r = modelInputs.r;
q = modelInputs.q;

rnCHF = modelInputs.rnCHF;
rnCHF = @(u)rnCHF(u).*exp(-1i*(r-q)*u*dt);  % undo the drift

%%%%%%%%%%%%%%%%%%%%%%%%%%
F_0 = gmdb_params.F_0;

g = gmdb_params.g;
c = gmdb_params.c;  % Cap on growth 
c_S = c;  % Cap on growth for surrender
gamma = gmdb_params.gamma;

alpha_fee = gmdb_params.alpha_fee;
death_prob = gmdb_params.death_prob_cond;

%%% NOTE: assume constant fee and r,q
lambda_m = log(1-alpha_fee) + (r-q)*dt;

%%%%%%%%%%%%%%%%%%%%%%%%

K = N/2;
dx = 2*alph/(N-1); a = 1/dx;

nnot = K/2;  % TODO: make sure 0 is on grid

%%%% Populate Beta coefficients for orthogonal projection
xmin = (1-K/2)*dx;

a2    = a^2;  
Cons2 = 24*a2*exp(-r*dt)/N;
zmin  = (1 - K)*dx;  %Kbar corresponds to zero

dw    = 2*pi*a/N;
grand = (dw: dw: (N-1)*dw);
grand = exp(-1i*zmin*grand).*rnCHF(grand).*(sin(grand/(2*a))./grand).^2./(2+cos(grand/a));
beta  = Cons2*real(fft([1/(24*a2) grand]));   %%%%  NOTE: all toep matrices incorporate exp(-r*dt)

toepM = [beta(K:-1:1)'; 0 ; beta(2*K-1:-1:K +1)'];
toepM = fft(toepM);

%%%% Initial terminal payoff coefficients (recursion proceeds backwards in time)
grid = xmin + dx*(0:K-1);
expGrid = exp(grid);

% Initialize Value function at maturity
Vals = F_0 * max(exp(g*T), min(exp(c*T), expGrid));

Thet = zeros(K,1);

cusum_right = cumsum(beta(2*K:-1:K +1))';
cusum_left = [ fliplr(cumsum(beta(1:1:K-1)))';0];

for m=M-1:-1:0

    % Define the value function at next time
    Vfunc = @(x)spline(grid, Vals, x);
    
    % Define maturity/death benefit function
    Mfunc = @(x) F_0 * max(exp(g*(m+1)*dt), min(exp(c*(m+1)*dt), exp(x)));
    
    % Combine the two functions, weight by death prob, and shift by lambda
    pw = death_prob(m+1);
    H = @(x) pw*Mfunc(lambda_m + x) + (1-pw)*Vfunc(lambda_m + x);
    
    Vals = H(grid);
    
    Thet(2:K -1) = (Vals(1:K-2)+10*Vals(2:K-1)+Vals(3:K))/12;
    Thet(1)      = (13*Vals(1)+15*Vals(2)-5*Vals(3)+Vals(4))/48;
    Thet(K)      = 2*(13*Vals(K)+15*Vals(K-1)-5*Vals(K-2)+Vals(K-3))/48;
   
    p = ifft(toepM.*fft([Thet(1:K);zeros(K,1)]));

    % update values with augmentation on the boundaries
    Vals(1:K) = p(1:K) + cusum_left*Vals(1) + cusum_right*Vals(end); 

    if m > 0 && gamma < 1  %%% Cant surrender in first period
        surrender = (1-gamma)*F_0*min(exp(c_S*m*dt), expGrid);
        Vals = max(Vals, surrender);
    end
end

price = Vals(nnot);

end

