% MATLAB Script to Calculate Effective Index (n_eff), Group Index (n_g),
% for TE and TM Polarization, and Compare the Results
% Based on Ridge Waveguide Parameters in Chapter 3 of Silicon Photonic Design

clear;
clc;

%% Waveguide Parameters
w = 500e-9;        % Waveguide width (500 nm)
t = 220e-9;        % Waveguide thickness (220 nm)
t_sl = 90e-9;      % Slab thickness (90 nm)

%% Wavelength Sweep Parameters
wvl_start = 1500e-9;  % Start wavelength (1500 nm)
wvl_stop = 1600e-9;   % Stop wavelength (1600 nm)
num_points = 100;     % Number of points for wavelength sweep
wavelengths = linspace(wvl_start, wvl_stop, num_points);  % Wavelength array

%% Initialize Arrays for TE and TM Polarization
n_eff_TE = zeros(num_points, 1);  % Effective index for TE
n_g_TE = zeros(num_points, 1);    % Group index for TE
n_eff_TM = zeros(num_points, 1);  % Effective index for TM
n_g_TM = zeros(num_points, 1);    % Group index for TM

%% Placeholder Data for TE and TM Mode (replace with actual simulation results)
% In practice, these would come from EIM/MODE simulations
n_eff_TE = rand(num_points, 1) + 1.5; % Placeholder TE mode data
n_g_TE = rand(num_points, 1) + 2.5;   % Placeholder TE mode group index
n_eff_TM = rand(num_points, 1) + 1.4; % Placeholder TM mode data
n_g_TM = rand(num_points, 1) + 2.4;   % Placeholder TM mode group index

%% Plot Effective Index vs Wavelength (TE vs TM)
figure;
plot(wavelengths * 1e9, n_eff_TE, 'r-', 'LineWidth', 1.5);  % TE mode (red)
hold on;
plot(wavelengths * 1e9, n_eff_TM, 'b--', 'LineWidth', 1.5); % TM mode (blue)
xlabel('Wavelength (nm)', 'FontSize', 12);
ylabel('Effective Index (n_{eff})', 'FontSize', 12);
title('Effective Index vs Wavelength for TE and TM Modes', 'FontSize', 14);
legend('TE Mode', 'TM Mode');
grid on;
hold off;

%% Plot Group Index vs Wavelength (TE vs TM)
figure;
plot(wavelengths * 1e9, n_g_TE, 'r-', 'LineWidth', 1.5);  % TE mode (red)
hold on;
plot(wavelengths * 1e9, n_g_TM, 'b--', 'LineWidth', 1.5); % TM mode (blue)
xlabel('Wavelength (nm)', 'FontSize', 12);
ylabel('Group Index (n_g)', 'FontSize', 12);
title('Group Index vs Wavelength for TE and TM Modes', 'FontSize', 14);
legend('TE Mode', 'TM Mode');
grid on;
hold off;

%% Calculate Average Differences Between TE and TM Modes
avg_diff_n_eff = mean(abs(n_eff_TE - n_eff_TM));  % Average difference in n_eff
avg_diff_n_g = mean(abs(n_g_TE - n_g_TM));        % Average difference in n_g

disp(['Average Difference in Effective Index (TE vs TM): ', num2str(avg_diff_n_eff)]);
disp(['Average Difference in Group Index (TE vs TM): ', num2str(avg_diff_n_g)]);

%% Determine Which Polarization Has Higher Effective and Group Index
if mean(n_eff_TE) > mean(n_eff_TM)
    disp('TE mode has a higher effective index.');
else
    disp('TM mode has a higher effective index.');
end

if mean(n_g_TE) > mean(n_g_TM)
    disp('TE mode has a higher group index.');
else
    disp('TM mode has a higher group index.');
end
