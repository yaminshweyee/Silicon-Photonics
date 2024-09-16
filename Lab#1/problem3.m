% MATLAB script to calculate effective index (n_eff), group index (n_g),
% and plot these for a silicon waveguide for TE and TM polarizations.

clear;
clc;

% Constants and waveguide parameters
w = 500e-9;        % Waveguide width in meters (500 nm)
t = 220e-9;        % Waveguide thickness in meters (220 nm)
t_sl = 90e-9;      % Slab thickness in meters (90 nm)

wvl_start = 1500e-9; % Starting wavelength in meters (1500 nm)
wvl_stop = 1600e-9;  % Stopping wavelength in meters (1600 nm)
num_points = 100;    % Number of points for wavelength sweep
wavelengths = linspace(wvl_start, wvl_stop, num_points); % Wavelength array

% Initialize arrays for storing results
n_eff_array_TE = zeros(num_points, 1);  % Effective index for TE
n_g_array_TE = zeros(num_points, 1);    % Group index for TE
n_eff_array_TM = zeros(num_points, 1);  % Effective index for TM
n_g_array_TM = zeros(num_points, 1);    % Group index for TM

% Assuming TE and TM results are loaded from files or obtained
% These variables should be populated with actual data
% For demonstration, using placeholder data:
n_eff_array_TE = rand(num_points, 1) + 1.5;  % Placeholder data for TE
n_g_array_TE = rand(num_points, 1) + 2.5;    % Placeholder data for TE
n_eff_array_TM = rand(num_points, 1) + 1.5;  % Placeholder data for TM
n_g_array_TM = rand(num_points, 1) + 2.5;    % Placeholder data for TM

% Plot Effective Index vs Wavelength for TE and TM
figure;
plot(wavelengths * 1e9, n_eff_array_TE, 'b-', 'LineWidth', 1.5); % TE in blue
hold on;
plot(wavelengths * 1e9, n_eff_array_TM, 'r--', 'LineWidth', 1.5); % TM in red
xlabel('Wavelength (nm)');
ylabel('Effective Index (n_{eff})');
title('Effective Index vs Wavelength');
legend('TE', 'TM');
grid on;
hold off;

% Plot Group Index vs Wavelength for TE and TM
figure;
plot(wavelengths * 1e9, n_g_array_TE, 'b-', 'LineWidth', 1.5); % TE in blue
hold on;
plot(wavelengths * 1e9, n_g_array_TM, 'r--', 'LineWidth', 1.5); % TM in red
xlabel('Wavelength (nm)');
ylabel('Group Index (n_g)');
title('Group Index vs Wavelength');
legend('TE', 'TM');
grid on;
hold off;

% Calculate and display average differences
avg_diff_n_eff = mean(abs(n_eff_array_TE - n_eff_array_TM));
avg_diff_n_g = mean(abs(n_g_array_TE - n_g_array_TM));

disp(['Average Difference in Effective Index (TE vs TM): ', num2str(avg_diff_n_eff)]);
disp(['Average Difference in Group Index (TE vs TM): ', num2str(avg_diff_n_g)]);

% Determine which polarization has a higher effective and group index
if mean(n_eff_array_TE) > mean(n_eff_array_TM)
    disp('TE polarization has a higher average effective index.');
else
    disp('TM polarization has a higher average effective index.');
end

if mean(n_g_array_TE) > mean(n_g_array_TM)
    disp('TE polarization has a higher average group index.');
else
    disp('TM polarization has a higher average group index.');
end
