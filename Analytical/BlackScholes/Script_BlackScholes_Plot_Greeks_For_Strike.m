%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descritpion: Script to plot Delta/Gamma of European option under Black Scholes Merton Model
% Author:      Justin Kirkby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------
%  Contract/Market Params
% ---------------------
call = 1;    %For call use 1 (else, its a put)
K    = 100;  %Initial price
r    = .10;  %Interest rate
q    = .00;  %dividend yield
T    = 1;    %Time (in years)

moneyness = [0.2:.01:1.8];
Svec = K * moneyness;   % strikes to price

sigma = 0.15;   % volatility

%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Prices (as function of moneyness)
h = figure();
subplot(5,1,1)
values = BSM_Greeks(0, S_0, sigma, r, q, T, Kvec, call);
plot(moneyness, values); grid on;
ylabel('price'); %xlabel('$K/S_0$', 'interpreter', 'latex')


% Plot Deltas (as function of moneyness)
subplot(5,1,2)
values = BSM_Greeks(1, S_0, sigma, r, q, T, Kvec, call);
plot(moneyness, values); grid on;
ylabel('delta'); %xlabel('$K/S_0$', 'interpreter', 'latex')

% Plot Gammas (as function of moneyness)
subplot(5,1,3)
values = BSM_Greeks(2, S_0, sigma, r, q, T, Kvec, call);
plot(moneyness, values); grid on;
ylabel('gamma');  % xlabel('$K/S_0$', 'interpreter', 'latex')


% Plot Thetas (as function of moneyness)
subplot(5,1,4)
values = BSM_Greeks(3, S_0, sigma, r, q, T, Kvec, call);
plot(moneyness, values); grid on;
ylabel('theta'); % xlabel('$K/S_0$', 'interpreter', 'latex')


% Plot Vegas (as function of moneyness)
subplot(5,1,5)
values = BSM_Greeks(4, S_0, sigma, r, q, T, Kvec, call);
plot(moneyness, values); grid on;
ylabel('vegas'); xlabel('$K/S_0$', 'interpreter', 'latex')


