%% Power-Law Differencing Scheme
clear; clc; close all;

% Parameters
L = 1.0;                    % Length of the domain (m)
phi0 = 100;                 % Boundary value at x = 0
phiL = 20;                  % Boundary value at x = L
rho = 1.0;                  % Fluid density (kg/m^3)
u = 1.00;                    % Velocity (m/s)
Gamma = 0.1;                % Diffusion coefficient
N = 11;                     % Number of grid points (increase for better resolution)

% Derived parameters
dx = L / (N - 1);           % Grid spacing
x = linspace(0, L, N);      % Spatial grid

% Initialize matrix A and vector b
A = zeros(N, N);            % Coefficient matrix
b = zeros(N, 1);            % RHS vector

% Build matrix A and RHS vector b
for i = 1:N
    if i == 1
        % Left boundary condition (Dirichlet)
        A(i, i) = 1;        % Set diagonal to 1
        b(i) = phi0;        % RHS equals boundary value
    elseif i == N
        % Right boundary condition (Dirichlet)
        A(i, i) = 1;        % Set diagonal to 1
        b(i) = phiL;        % RHS equals boundary value
    else
        % Internal nodes
        Local_Pe = rho * u * dx / Gamma; % Local Péclet number
        f_Pe = max((1 - 0.1 * abs(Local_Pe))^5, 0); % Power-law weighting function

        alpha_W = (Gamma / dx) * f_Pe + max(rho * u, 0);   % West coefficient
        alpha_E = (Gamma / dx) * f_Pe + max(-rho * u, 0);  % East coefficient
        alpha_P = alpha_W + alpha_E;                    % Central coefficient

        A(i, i-1) = -alpha_W;    % West coefficient
        A(i, i) = alpha_P;       % Central coefficient
        A(i, i+1) = -alpha_E;    % East coefficient
        b(i) = 0;                % No source term
    end
end

% Solve the system using matrix inversion
phi = A \ b;

% Analytical solution 
Pe = rho * u * L / Gamma; % Péclet number

phi_analytical = (phi0 + (phiL - phi0) * (exp(Pe * x / L) - 1) ./ (exp(Pe) - 1))';

% Compute percentage error
error = abs((phi - phi_analytical) ./ (phi_analytical)) .* 100;
Error = sum(error);

% Plot results
figure;
plot(x, phi, 'bo--','LineWidth', 0.5 , 'MarkerFaceColor', 'b', 'DisplayName', 'Numerical (PLDS)');
hold on;
plot(x, phi_analytical, 'r-', 'DisplayName', 'Analytical Solution');
xlabel('\it x (m)','Fontsize', 15, 'Interpreter','latex');
ylabel('\phi');
set(gca, 'YDir', 'reverse'); % Invert the y-axis
ylim([20 120]); % Set y-axis limits
xlim([0 1]); % Ensure x-axis ranges from 0 to 1
legend('Location', 'SouthWest'); % Adjust legend position
grid on;
title('Analytical vs Numerical solution using PLDS','Fontsize', 15 ,'Interpreter','latex');
subtitle(['\rho = ' num2str(rho) ' kg/m^{3}', '    \Gamma_{\phi} = ' num2str(Gamma) 'm^{2}/s', '    \it Local Pe = ' num2str(Local_Pe),  '    \it Global Pe = ' num2str(Pe),  '    \it \color{blue} u = ' num2str(u) 'm/s'],'Fontsize', 15)
text(-0.07, 100, '\phi_{0} = 100','Fontsize', 10);
text(-0.06, 21, '\phi_{L} = 20','Fontsize', 10)
errortext = ['Error =' num2str(Error) '%'];
text(0.5, 50, errortext,"FontSize", 15);
hold off;

% Display the error
fprintf('Mean Percentage Error: %.2f%%\n', Error);
