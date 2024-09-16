% MATLAB script to calculate effective index (n_eff), group index (n_g),
% and plot these for a silicon waveguide based on Chapter 3 of Silicon Photonic Design.

clear;
clc;

% Constants and waveguide parameters
w = 500e-9;        % waveguide width in meters
t = 220e-9;        % waveguide height (thickness) in meters
wvl_start = 1500e-9; % start wavelength in meters
wvl_stop = 1600e-9;  % stop wavelength in meters
num_points = 100;    % number of points for wavelength sweep
wavelengths = linspace(wvl_start, wvl_stop, num_points); % wavelength range

% Dispersive refractive indices for Silicon and SiO2
n_si = @(lambda) 3.48 - 0.0005 * (lambda - 1550e-9) * 1e6; % Silicon index approx.
n_sio2 = @(lambda) 1.444 + 0.0001 * (lambda - 1550e-9) * 1e6; % SiO2 index approx.

% Approximation for the Effective Index Method (EIM) for a silicon waveguide
effective_index_EIM = @(lambda) n_si(lambda) - 0.05;  % Approximate value for EIM

% Numerical MODE solver approximation (can be replaced with actual simulation)
effective_index_MODE = @(lambda) n_si(lambda) - 0.07;  % More accurate with MODE solver

% Function to calculate group index from effective index and wavelength
group_index = @(n_eff, lambda) n_eff - lambda .* gradient(n_eff, lambda);

% Calculate effective indices for both methods
n_eff_EIM = effective_index_EIM(wavelengths);
n_eff_MODE = effective_index_MODE(wavelengths);

% Calculate group indices for both methods
n_g_EIM = group_index(n_eff_EIM, wavelengths);
n_g_MODE = group_index(n_eff_MODE, wavelengths);

% Dispersion Calculation: D(lambda) = -(lambda/c) * d²n_eff/dlambda²
c = 3e8; % Speed of light in vacuum (m/s)
d2n_eff_dlambda2_EIM = gradient(gradient(n_eff_EIM, wavelengths), wavelengths);
d2n_eff_dlambda2_MODE = gradient(gradient(n_eff_MODE, wavelengths), wavelengths);

dispersion_EIM = -(wavelengths / c) .* d2n_eff_dlambda2_EIM;
dispersion_MODE = -(wavelengths / c) .* d2n_eff_dlambda2_MODE;

% Plot Effective Index vs Wavelength
figure;
plot(wavelengths * 1e9, n_eff_EIM, '--', 'LineWidth', 1.5);
hold on;
plot(wavelengths * 1e9, n_eff_MODE, '-', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Effective Index (n_{eff})');
title('Effective Index vs Wavelength');
legend('EIM Effective Index', 'MODE Effective Index');
grid on;

% Plot Group Index vs Wavelength
figure;
plot(wavelengths * 1e9, n_g_EIM, '--', 'LineWidth', 1.5);
hold on;
plot(wavelengths * 1e9, n_g_MODE, '-', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Group Index (n_g)');
title('Group Index vs Wavelength');
legend('EIM Group Index', 'MODE Group Index');
grid on;

% Plot Dispersion vs Wavelength
figure;
plot(wavelengths * 1e9, dispersion_EIM, '--', 'LineWidth', 1.5);
hold on;
plot(wavelengths * 1e9, dispersion_MODE, '-', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Dispersion (s/m^2)');
title('Dispersion vs Wavelength');
legend('EIM Dispersion', 'MODE Dispersion');
grid on;

% Calculate average differences for effective index and group index
avg_diff_n_eff = mean(abs(n_eff_EIM - n_eff_MODE));
avg_diff_n_g = mean(abs(n_g_EIM - n_g_MODE));

% Display the average differences
fprintf('Average Difference in Effective Index: %.4f\n', avg_diff_n_eff);
fprintf('Average Difference in Group Index: %.4f\n', avg_diff_n_g);
