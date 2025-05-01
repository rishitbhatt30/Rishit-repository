%% Upwind Differencing Scheme: Error Analysis for Multiple u Values
clear; clc; close all;

% Parameters
L = 1.0;                    % Length of the domain (m)
phi0 = 100;                 % Boundary value at x = 0
phiL = 20;                  % Boundary value at x = L
rho = 1.0;                  % Fluid density (kg/m^3)
Gamma = 0.1;                % Diffusion coefficient
u_values = [1, 3, 5];       % Different velocities to analyze

% Colors for plotting
colors = {'r-', 'g-', 'b-'};
legend_labels = {};

% Initialize figure for plotting
figure;
hold on;

for k = 1:length(u_values)
    u = u_values(k);        % Current velocity value
    dx_values = [];         % Array to store dx values
    Error_values = [];      % Array to store Error values

    for N = 7:31            % Number of grid points
        % Derived parameters
        dx = L / (N - 1);   % Grid spacing
        x = linspace(0, L, N);  % Spatial grid

        % Initialize matrix A and vector b
        A = zeros(N, N);    % Coefficient matrix
        b = zeros(N, 1);    % RHS vector

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
                alpha_W = Gamma / dx + max(rho * u, 0);   % West coefficient
                alpha_E = Gamma / dx + max(-rho * u, 0);  % East coefficient
                alpha_P = alpha_W + alpha_E;             % Central coefficient

                A(i, i-1) = -alpha_W;    % West coefficient
                A(i, i) = alpha_P;       % Central coefficient
                A(i, i+1) = -alpha_E;    % East coefficient
                b(i) = 0;                % No source term
            end
        end

        % Solve the system using matrix inversion
        phi = A \ b;

        % Analytical solution 
        Pe = rho * u * L / Gamma; % Global Péclet number
        phi_analytical = (phi0 + (phiL - phi0) * (exp(Pe * x / L) - 1) ./ (exp(Pe) - 1))';

        % Compute percentage error
        error = abs((phi - phi_analytical) ./ (phi_analytical)) .* 100;
        total_error = sum(error); % Summed error

        % Store dx and total_error
        dx_values = [dx_values; dx];
        Error_values = [Error_values; total_error];
    end

    % Plot dx vs Error for the current u
    loglog(dx_values, Error_values, colors{k}, 'LineWidth', 1.5);
    legend_labels{k} = sprintf('u = %d', u); % Add legend label
end

% Graph formatting
xlabel('Grid Spacing (\Delta x) / m');
ylabel('Total Error (%)');
grid on;
title('Error vs Grid Spacing for UDS at Different u Values');
legend(legend_labels, 'Location', 'SouthWest');
hold off;
