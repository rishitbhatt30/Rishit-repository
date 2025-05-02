
%%

% Lubricant Viscosity Calculation with Row Number Selection
% Equation: log(log(η + 0.7)) = A - B * log(T in Kelvin)

clc; clear; close all;

% === Load Lubricant Data from Excel ===
filename = 'Lubricant_Viscosity_Table.xlsx'; 
sheet = 1;
data = readtable(filename, 'Sheet', sheet);

% === Display Available Lubricants with Row Numbers ===
disp('Available Lubricants:');
disp(table((1:height(data))', data.Lubricant, 'VariableNames', {'Index', 'Lubricant'}));

% === User Input: Select Lubricant by Row Number ===
row_idx = input('Enter the row number of the lubricant: ');

if row_idx < 1 || row_idx > height(data)
    error('Invalid row number. Please enter a valid row.');
 
end

% === Extract Viscosity Values from Table ===
ISO_VG = data.ISO_VG(row_idx);      % Column 2
eta40 = data{row_idx, 4};           % Column 4 (Viscosity at 40°C)
eta100 = data{row_idx, 5};          % Column 5 (Viscosity at 100°C)
k = data{row_idx, 6};               % Column 6 (Pressure-viscosity coefficient k)
s = data{row_idx, 7};               % Column 7 (Exponent s)

% === User Input: Operating Temperature ===
T_op_C = input('Enter operating temperature (°C): ');
T_op = T_op_C + 273.15;  % Convert to Kelvin

% === Solve for A and B (Using Temperature in Kelvin) ===
T40 = 40 + 273.15;  % Convert 40°C to Kelvin
T100 = 100 + 273.15; % Convert 100°C to Kelvin

loglog_eta40 = log(log(eta40 + 0.7)); % log(log(η + 0.7)) at 40°C
loglog_eta100 = log(log(eta100 + 0.7)); % log(log(η + 0.7)) at 100°C
logT40 = log(T40);
logT100 = log(T100);

% Solve for B
B = (loglog_eta40 - loglog_eta100) / (logT40 - logT100);

% Solve for A
A = loglog_eta40 + B * logT40;

% === Calculate Viscosity at Given Temperature ===
loglog_etaTop = A - B * log(T_op);
eta_Top = exp(exp(loglog_etaTop)) - 0.7;  % Solve for η
eta_Top_Pas = eta_Top* 850;  % Convert cSt to Pa·s (Assuming density = 850 kg/m³)

% === Calculate Pressure-Viscosity Coefficient ===
alpha_1 = k * (eta_Top_Pas)^s;  % in Pa^-1

% === Display Results ===
fprintf('\nResults:\n');
fprintf('Selected Lubricant: %s\n', data.Lubricant{row_idx});
fprintf('ISO VG: %d\n', ISO_VG);
fprintf('Viscosity at %d°C (%.2f K): %.2f mm²/s\n', T_op_C, T_op, eta_Top);
fprintf('Pressure-Viscosity Coefficient (α₁): %.2e Pa⁻¹\n', alpha_1);

% === Check if viscosity meets AGMA Guidelines ===
if eta_Top >= 15
    fprintf('✔ The lubricant is suitable for heavy loads.\n');
elseif eta_Top >= 8
    fprintf('✔ The lubricant is suitable for moderate loads.\n');
else
    fprintf('⚠ Warning: The viscosity is too low for reliable lubrication!\n');
end

% === Plot Viscosity-Temperature Curve ===
T_range_C = linspace(30, 120, 100); % Temperature range from 30°C to 120°C
T_range_K = T_range_C + 273.15;  % Convert to Kelvin
loglog_eta_range = A - B * log(T_range_K);
eta_range = exp(exp(loglog_eta_range)) - 0.7;

figure;
semilogy(T_range_C, eta_range, 'b', 'LineWidth', 2); hold on;
scatter(T_op_C, eta_Top, 100, 'r', 'filled'); % Mark operating temperature
xlabel('Temperature (°C)');
ylabel('Viscosity (mm²/s)');
title(sprintf('Viscosity-Temperature Curve (%s, ISO VG %d)', data.Lubricant{row_idx}, ISO_VG));
grid on;
legend('Viscosity Curve', 'Operating Point');


