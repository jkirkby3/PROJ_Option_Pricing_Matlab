%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to calc prices and Greeks European options under Black Scholes Merton Model
% Author:      Justin Kirkby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------
%  Contract/Market Params
% ---------------------
call = 1;    %For call use 1 (else, its a put)
S_0  = 100;  %Initial price
r    = .00;  %Interest rate
q    = .00;  %dividend yield
T    = 1;    %Time (in years)
Kvec = S_0 * [0.2:.01:1.8];   % strikes to price

sigma = 0.2;   % volatility

% ---------------------
% Select what to calculate
% G = Greek: 
%         0:Price
%         1:Delta, 2:Gamma, 3:Theta
%         4:Vega,  5:Rho,   6:Vanna
%         7:Vomma
%  
% ---------------------
G = 0;   
% ---------------------
         
values = BSM_Greeks(G, S_0, sigma, r, q, T, Kvec, call);


% Plot
plot(Kvec, values)
ylabel('price (or greek)')
xlabel('strike')
grid on;

