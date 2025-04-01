% General share info
S_0 = [100; 95; 50];  % Initial stock prices
sigma = [0.15; 0.2; 0.3];  % Volatilities
corr_mat = [1, 0.2, 0.4; 0.2, 1, 0.8; 0.4, 0.8, 1];  % Correlation matrix
L = chol(corr_mat, 'lower');  % Cholesky decomposition
r = 0.1;  % Risk-free rate
T = 1;  % Time horizon in years
t_simulations = 10000;  % Number of simulations
alpha = 0.05;  % Significance level for VaR calculation

% Current portfolio value
portval_current = sum(S_0);

% Creating 10000 simulations for future portfolio values
Z = L * normrnd(0, 1, [3, t_simulations]);  % Generate random samples
portval_future = sum(terminal_shareprice(S_0, r, sigma, Z, T), 1);  % Calculate future portfolio values

% Calculating portfolio returns
port_return = (portval_future - portval_current) / portval_current;

% Sorting the Returns
port_return = sort(port_return);

% Determining VaR (Value at Risk)
mVar_estimate = -port_return(floor(alpha * t_simulations) + 1);  % VaR at alpha confidence level

% Display the result
disp(['Value at Risk (VaR) at alpha = ', num2str(alpha), ' is: ', num2str(mVar_estimate)]);

%% Function to calculate terminal share price
function terminal_price = terminal_shareprice(S_0, risk_free_rate, sigma, Z, T)
    % This function generates the terminal share price given some random normal values Z
    terminal_price = S_0 .* exp((risk_free_rate - sigma.^2 / 2) * T + sigma .* sqrt(T) .* Z);
end
