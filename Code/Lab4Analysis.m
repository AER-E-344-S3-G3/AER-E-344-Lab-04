% AER E 344 Spring 2024 Lab 04 Analysis
% Section 3 Group 3
clear, clc, close all;

%% Constants
u = symunit;
figure_dir = "../Figures/";
data_zip = "./data.zip";
K = 1.1; % []
% Viscosity calculated at 21.2°C
% https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601
% .html
mu_air = 18.18e-6; % [Pa*s]
rho_air = 1.225; % [kg/m^3]
p_A_idx = 19;
p_E_idx = 20;
delta_theta = deg2rad(20);
angles = 180 : -20 : -180; % [°]
angles_shortened = 180 : -20 : -160; % [°]
D = double(separateUnits(unitConvert(0.7 * u.in, u.m))); % [m]

%% Import Calibration Data
unzip(data_zip);
freq = 5 : 5 : 35; % [Hz]
freq_str = strings([1 length(freq)]);
for i = 1 : length(freq_str)
    freq_str(i) = freq(i) + " Hz";
end
%{
Pressure Matrix Structure:
Each row is a different frequency corresponding to the values in the freq
array. Each column is a different pressure tap. For example, to get the 6th
pressure tap for 10 Hz, you would call pressure_matrix(2, 6).
%}
p = zeros([length(freq) 20]); % [Pa]
% defines the columns we actually want to parse
column_mask = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 1 1 1 1];

for i = 1 : length(freq)
    data = readmatrix(sprintf("%02dHz.csv", freq(i)));
    counter = 1;
    for j = 1 : length(column_mask)
        % Only find the mean of certain columns
        if column_mask(j)
            p(i, counter) = mean(data(:, j)); % [Pa]
            counter = counter + 1;
        end
    end
end
clear data;

%% Calculations
C_p = zeros([length(freq) 18]);
C_d = zeros([1 length(freq)]); % []
Re = zeros([1 length(freq)]); % []
for i = 1 : length(freq)
    C_p(i, :) = (p(i, 1:18) - p(i, p_E_idx)) ...
        ./ (K .* (p(i, p_A_idx) - p(i, p_E_idx))); % []
    C_d(i) = -delta_theta ./ 2 ...
        .* sum(C_p(i, :) .* cosd(angles_shortened)); % []
    V = sqrt(K * (p(i, p_A_idx) - p(i, p_E_idx)) * 2 / rho_air);
    Re(i) = rho_air * V * 2 * D / mu_air; % []
end

% Make the data symmetrical
C_d_trap = zeros([1 length(freq)]); % []
for i = 1 : length(freq)
    C_p(i, 19) = C_p(i, 1); % []
    C_d_trap(i) = -delta_theta / 2 * trapz(C_p(i, :) .* cosd(angles)); % []
end

theta_theoretical = 0 : 0.01 : 360; % [°]
C_p_theoretical = 1 - 4 .* (sind(theta_theoretical).^2); % []

%% Plot
% Plot the individual C_p graphs
figure_cnt = 1;
for i = 1 : length(freq)
    figure(figure_cnt);
    plot(0 : 20 : 360, C_p(i, :), "-x", "MarkerSize", 8, "LineWidth", 2);
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = ...
        sprintf("C_P Distribution of a Cylinder at %d Hz", freq(i));
    title(title_str);
    xlabel("Angle [°]");
    ylabel("C_P [ ]");
    legend("Coefficient of Pressure, C_P");
    xlim([0 360]);
    grid on;
    figure_cnt = figure_cnt + 1;
    saveas(gcf, figure_dir + title_str + ".svg");
end

% Plot a cumulative C_p graph
figure(figure_cnt);
p = plot(theta_theoretical, C_p_theoretical, "--", "LineWidth", 1.5);
hold on
for i = 1 : length(freq)
    plot(0 : 20 : 360, C_p(i, :))
end
hold off;
fontname("Times New Roman");
fontsize(12, "points");
title_str = "C_P Distribution of a Cylinder";
title(title_str);
xlabel("Angle [°]");
ylabel("C_P [ ]");
legend(["Theoretical Solution" freq_str]);
xlim([0 360]);
uistack(p, "bottom");
grid on;
figure_cnt = figure_cnt + 1;
saveas(gcf, figure_dir + title_str + ".svg");

% Plot a cumulative C_p graph without 5-15 Hz
figure(figure_cnt);
p = plot(theta_theoretical, C_p_theoretical, "--", "LineWidth", 1.5);
hold on;
for i = 4 : length(freq)
    plot(0 : 20 : 360, C_p(i, :));
end
hold off;
fontname("Times New Roman");
fontsize(12, "points");
title_str = "C_P Distribution of a Cylinder (Excluding 5–15 Hz)";
title(title_str);
xlabel("Angle [°]");
ylabel("C_P [ ]");
legend(["Theoretical Solution" freq_str(4 : end)]);
xlim([0 360]);
uistack(p, "bottom");
grid on;
figure_cnt = figure_cnt + 1;
saveas(gcf, figure_dir + title_str + ".svg");

% Plot the C_d graph
figure(figure_cnt);
plot(freq, C_d, ".-", "MarkerSize", 20, "LineWidth", 3);
hold on;
plot(freq, C_d_trap, ".--", "MarkerSize", 20, "LineWidth", 2);
hold off;
fontname("Times New Roman");
fontsize(12, "points");
title_str = "C_d of a Cylinder";
title(title_str);
xlabel("\omega_{motor} [Hz]");
ylabel("C_d [ ]");
legend("Riemann Sum", "Trapezoidal Integration");
grid on;
figure_cnt = figure_cnt + 1;
saveas(gcf, figure_dir + title_str + ".svg");

% Plot the Reynold's graph
figure(figure_cnt);
plot(freq, Re, ".-", "MarkerSize", 20, "LineWidth", 2);
fontname("Times New Roman");
fontsize(12, "points");
title_str = "Reynolds Number for a Cylinder";
title(title_str);
xlabel("\omega_{motor} [Hz]");
ylabel("Re [ ]");
grid on;
figure_cnt = figure_cnt + 1;
saveas(gcf, figure_dir + title_str + ".svg");

% Plot C_d vs Reynold's
figure(figure_cnt);
plot(Re, C_d, ".-", "MarkerSize", 20, "LineWidth", 2);
fontname("Times New Roman");
fontsize(12, "points");
title_str = "C_d vs. Reynolds Number for a Cylinder";
title(title_str);
xlabel("Re []");
ylabel("C_d [ ]");
grid on;
figure_cnt = figure_cnt + 1;
saveas(gcf, figure_dir + title_str + ".svg");

delete *.csv