% Define waveguide parameters
w = 500e-9;       % Width of the waveguide in meters (500 nm)
t = 220e-9;       % Thickness of the waveguide in meters (220 nm)
t_sl = 90e-9;     % Slab thickness in meters (90 nm)

% Wavelength sweep parameters
wvl_start = 1500e-9;  % Starting wavelength in meters (1500 nm)
wvl_stop = 1600e-9;   % Stopping wavelength in meters (1600 nm)
num_points = 100;     % Number of points for the wavelength sweep
wavelengths = linspace(wvl_start, wvl_stop, num_points);  % Wavelength array

% Initialize arrays for storing results
n_eff_array_EIM = zeros(num_points, 1);  % Effective index for EIM
n_g_array_EIM = zeros(num_points, 1);    % Group index for EIM
n_eff_array_MODE = zeros(num_points, 1); % Effective index for MODE
n_g_array_MODE = zeros(num_points, 1);   % Group index for MODE

% Assuming EIM and MODE results are loaded from files or obtained
% These variables should be populated with actual data
% For demonstration, using placeholder data:
n_eff_array_EIM = rand(num_points, 1) + 1.5;  % Placeholder data
n_g_array_EIM = rand(num_points, 1) + 2.5;    % Placeholder data
n_eff_array_MODE = rand(num_points, 1) + 1.5; % Placeholder data
n_g_array_MODE = rand(num_points, 1) + 2.5;   % Placeholder data

% Plot Effective Index vs Wavelength
figure;
plot(wavelengths * 1e9, n_eff_array_EIM, 'r-', 'LineWidth', 1.5); % EIM in red
hold on;
plot(wavelengths * 1e9, n_eff_array_MODE, 'b--', 'LineWidth', 1.5); % MODE in blue
xlabel('Wavelength (nm)');
ylabel('Effective Index (n_{eff})');
title('Effective Index vs Wavelength');
legend('EIM', 'MODE');
grid on;
hold off;

% Plot Group Index vs Wavelength
figure;
plot(wavelengths * 1e9, n_g_array_EIM, 'r-', 'LineWidth', 1.5); % EIM in red
hold on;
plot(wavelengths * 1e9, n_g_array_MODE, 'b--', 'LineWidth', 1.5); % MODE in blue
xlabel('Wavelength (nm)');
ylabel('Group Index (n_g)');
title('Group Index vs Wavelength');
legend('EIM', 'MODE');
grid on;
hold off;

% Calculate and display average differences
avg_diff_n_eff = mean(abs(n_eff_array_EIM - n_eff_array_MODE));
avg_diff_n_g = mean(abs(n_g_array_EIM - n_g_array_MODE));

disp(['Average Difference in Effective Index: ', num2str(avg_diff_n_eff)]);
disp(['Average Difference in Group Index: ', num2str(avg_diff_n_g)]);
