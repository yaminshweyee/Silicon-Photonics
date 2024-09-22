
% MATLAB script to sweep waveguide width and plot effective index (n_eff) and group index (n_g)

clear;
clc;

% Constants and parameters
wvl = 1550e-9; % Optical wavelength in meters (1550 nm)
width_start = 300e-9; % Starting width (300 nm)
width_stop = 1200e-9; % Ending width (1200 nm)
step_size = 5e-9; % Step size (5 nm)
num_points = (width_stop - width_start) / step_size + 1; % Number of points

% Initialize arrays to store results
waveguide_widths = linspace(width_start, width_stop, num_points); % Array of waveguide widths
n_eff_values = zeros(1, num_points); % Array to store effective index values
n_g_values = zeros(1, num_points); % Array to store group index values

% Loop through each width
for i = 1:num_points
    width = waveguide_widths(i); % Current waveguide width
    
    % You would call the MODE solver here for each width.
    % For this example, we simulate the results with a function.
    % Replace this with actual MODE solver calls.
    
    % Example of calculating effective and group indices (placeholder formulas):
    n_eff = 2.8 + 0.01 * randn(); % Random noise added to simulate slight variation in n_eff
    n_g = 3.8 + 0.01 * randn(); % Random noise added to simulate slight variation in n_g
    
    % Store the results
    n_eff_values(i) = n_eff;
    n_g_values(i) = n_g;
end

% Plot Effective Index vs Waveguide Width
figure;
plot(waveguide_widths * 1e9, n_eff_values, '-o', 'LineWidth', 1.5);
xlabel('Waveguide Width (nm)');
ylabel('Effective Index (n_{eff})');
title('Effective Index vs Waveguide Width');
grid on;

% Plot Group Index vs Waveguide Width
figure;
plot(waveguide_widths * 1e9, n_g_values, '-o', 'LineWidth', 1.5);
xlabel('Waveguide Width (nm)');
ylabel('Group Index (n_g)');
title('Group Index vs Waveguide Width');
grid on;

% Save results to file (optional)
save('waveguide_width_sweep_results.mat', 'waveguide_widths', 'n_eff_values', 'n_g_values');

% Display average results
fprintf('Average Effective Index: %.4f\n', mean(n_eff_values));
fprintf('Average Group Index: %.4f\n', mean(n_g_values));
