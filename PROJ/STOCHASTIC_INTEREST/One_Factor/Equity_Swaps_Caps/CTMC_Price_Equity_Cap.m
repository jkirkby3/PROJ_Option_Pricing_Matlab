function [ price ] = CTMC_Price_Equity_Cap( numeric_param, M, dt, model, modparam, contract_param )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About: Price Equity Cap/Floors under 1 Factor short rate model
% Models Supported: 1 Factor short rate model
% Returns: price of bond
% Author: Justin Lars Kirkby
%
% References:  (1) 
%
% ----------------------
% Contract Params 
% ----------------------
% T  : time remaining to maturity (T = 2 is two years, T = .5 is half a year)
%
% ----------------------
% Model Params 
% ----------------------
% model:
%        1 = HESTON:      Sigmar, r0, rho, eta, theta
%        2 = VASICEK:     Sigmar, r0, rho, eta, theta
%
% modparam: contains all necessary params for the specific model (see below during assingment which ones are needed)
%
% ----------------------
% Numerical Params 
% ----------------------
% numeric_parm: container of numerical params
%   m_0: number of states to approximate the Heston model with
%   gamma: var grid width parameter, grid is +/- gamma*stddev(variance process)
%   gridMethod: which type of var grid to use (typcially use 4)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_0 = numeric_param.m_0;
gridMethod = numeric_param.gridMethod;
gamma = numeric_param.gamma;
gridMultParam = numeric_param.gridMultParam;

N = numeric_param.N;

contract_type = contract_param.contract_type; % 1 = swap, 2 = cap/floor
K = contract_param.K;
F = contract_param.F;
C = contract_param.C;
c = contract_param.c;

rho = modparam.rho;
sigma = modparam.Sigma_s;

Tmax = M*dt;
t = Tmax/3;
[lx, r0, ux] = get_rate_grid_boundaries( model, modparam, t, gamma);

[ mu_func,  sig_func, zeta_func, varrho_func, nu, beta_, is_affine_varrho]  ...
    = get_short_rate_grid_diffusion_funcs( model,  modparam);

boundaryMethod = 1;
center = r0; %this is where grid clusters... we can experiment with other choices.. 

[Q, rs, c_index] = Q_Matrix(m_0, mu_func,sig_func,lx,ux,gridMethod,center, gridMultParam, boundaryMethod);
Qt = dt*Q;

Z = expm(Qt - dt*diag(rs));
gam_tilde = make_gamma_tilde(rs, dt, c);  % NOTE: this is 1 + r_til_{i,j} + c

experimental = 0;

if contract_type == 1 && experimental == 0
    % Do efficient approach
    P = expm(Qt);
    E = P - Z.*gam_tilde;  % "payoff" matrix

elseif contract_type == 2 || (contract_type == 1 && experimental == 1)
    if is_affine_varrho ~= 1 && rho ~= 0
        error("This model is not affine, only rho=0 is supported");
    end
    
    E = zeros(m_0, m_0);
    P = expm(Qt);
    
    sigrho2 = sigma^2*(1 - rho^2)*dt;
    sigrho = sqrt(sigrho2);
    cons1 = (nu + sigma^2/2)*dt;
    
    %%% PROJECTION

    m_l = rs(1)*dt;
    m_u = rs(m_0)*dt;
    dx = (m_u - m_l) / (N - 1);  a = 1/dx;
    
    Cons = 24*a^2/N;
    Cons2 = a^0.5*Cons;
    
    dxi    = 2*pi*a/N;
    xi     = dxi*(0:N-1)';
    
    EXP_A  = ones(m_0,m_0,N-1);    % THis will be matrix exponential of integrated variance
    for n = 1:N  
        EXP_A(:,:,n) = expm(Qt + (xi(n)*1i*dt)*diag(rs));  %This incorporates the dt   
    end
    
    % Determine the point masses, and remove from continuous portion of cdf
    pj = zeros(m_0,1);
    for k = 1:m_0
        cj = rs(k)*dt;
        pj(k) = 1 - exp(Qt(k,k));  % point mass, not that Qt(kk) = dt*Q(kk) along diagonal
        EXP_A(k,k,:) = ( EXP_A(k,k,:) - (1 - pj(k))*reshape(exp(1i*xi*cj), 1, 1, N) ) / pj(k);
    end
    %plot(1-pj)
    
    xmin = m_l;
    grand_xi = exp(-1i*xmin*xi).*(sin(xi/(2*a))./xi).^2./(2+cos(xi/a));
    grand_xi = grand_xi(2:N);
    
    ys = m_l + (dx/2)*(0:2*(N-1));
    
    tol_prob = 1e-12;
    
    
    for k = 1:m_0
        
        for l_=1:m_0
            if P(k,l_) < tol_prob
                continue;
            end
         
            zetakl = zeta_func(rs(k), rs(l_));
            %zetakl = zetakl - lamk;
            
            mukly = (zetakl - cons1) + ys*(1-beta_);       
            gamkl = gam_tilde(k,l_);
            
            % Define H(y)*exp(-y)
            if contract_type == 1
                H = exp(-ys).*(exp(mukly + sigrho2/2) - gamkl);
            else
                %%% Define the constants (dont depend on y)
                aklc = log(max(0, F + gamkl)) / sigrho;
                bklc = log(max(0, C + gamkl)) / sigrho;
                
                akly = aklc - mukly / sigrho;
                bkly = bklc - mukly / sigrho;
                
                nb = normcdf(bkly); na = normcdf(akly);
                
                H = exp(mukly + sigrho2/2) ...
                    .* (normcdf(bkly - sigrho) - normcdf(akly - sigrho)) ...
                    - gamkl * (nb - na) + F*na + C*(1-nb);
                H = exp(-ys).*H;
            end

            
            %%% Calculate Fourier Inversion to get PDF

            if abs(EXP_A(k,l_,1)) > 1e-12
                % NOTE: EXP_A(k,l_,1) is the probalility AFTER removing
                % point mass

                beta = (Cons2/EXP_A(k,l_,1)) ...
                            *real(fft([EXP_A(k,l_,1)/(24*a^2);
                            grand_xi.*reshape(EXP_A(k,l_,2:N), N-1, 1)]));

				% NOTE: 1 = sum(beta)*a^-0.5
                thetas = zeros(size(beta));
                thetas(1) = (H(1) + 2*H(2))/3;
                thetas(N) = (2*H(N-1) + H(N))/3;
                thetas(2:N-1) = (H(2:2:2*N-3) + H(3:2:2*N-2) + H(4:2:2*N-1))/3;
                E(k,l_) = a^-0.5*sum(thetas.*beta);
            end
  
            % Numerical Integral
            
            if k==l_
                ck = rs(k)*dt;
                
                muu = (zetakl - cons1) + ck * (1-beta_);
                if contract_type == 1 
                    Hck = exp(-ck)*(exp(muu + sigrho2/2) - gamkl);
                else
                    %%% Define the constants (dont depend on y)
                    akly = aklc - muu / sigrho;
                    bkly = bklc - muu / sigrho;
                    nb = normcdf(bkly); na = normcdf(akly);

                    H = exp(muu + sigrho2/2) ...
                        .* (normcdf(bkly - sigrho) - normcdf(akly - sigrho)) ...
                        - gamkl * (nb - na) + F*na + C*(1-nb);
                    Hck = exp(-ck)*H;
                end
                
                E(k,l_) =  pj(k)*E(k,l_) + (1-pj(k))*exp(-ck)*Hck;
            end
            
            E(k,l_) = P(k,l_)*E(k,l_);
        end
    end
end

% Compute the sum over integrals (indendent of E)
H = sum(E, 2);

%%% Calculate price
V = H;
for m=2:M
   V = H + Z*V; 
end

price = K*V(c_index);

end

function [gam_tilde] = make_gamma_tilde(rs, dt, c)
m_0 = length(rs);
gam_tilde = zeros(m_0, m_0);  % NOTE: this is 1 + r_til_{i,j} + c

% Todo: more efficient
for i = 1:m_0
    for j=1:m_0
        gam_tilde(i,j) = dt*(rs(i) + rs(j))/2 + 1 + c;
    end
end

end


