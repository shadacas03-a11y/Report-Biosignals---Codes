%% === Poincaré Plot Analysis (Full + Zoomed) ===
clc; clear; close all;

%% === Load ECG Data ===
filename = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';
data = readmatrix(filename);
fs = 100;                      % Sampling frequency [Hz]
ECG = data(:,6);               % Lead II
t = (0:length(ECG)-1)/fs;

%% === Filter ECG Signal ===
[b_hp,a_hp] = butter(2,0.5/(fs/2),'high');   % High-pass
[b_lp,a_lp] = butter(4,35/(fs/2),'low');    % Low-pass
ECG_filt = filtfilt(b_lp,a_lp,filtfilt(b_hp,a_hp,ECG))*1000;  % mV
ECG_smooth = smoothdata(ECG_filt,'gaussian',round(0.015*fs));

%% === Detect R Peaks ===
minPeakHeight = 0.25*max(ECG_smooth);
minPeakDistance = 0.5*fs;
[R_peaks, R_locs] = findpeaks(ECG_smooth,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);

%% === Compute RR Intervals and Remove Outliers ===
RR_intervals = diff(t(R_locs));
remove_outliers = @(x) x(x >= prctile(x,25)-1.5*(prctile(x,75)-prctile(x,25)) & ...
                          x <= prctile(x,75)+1.5*(prctile(x,75)-prctile(x,25)));
RR_clean = remove_outliers(RR_intervals);

%% === Poincaré Full Data ===
RR_n   = RR_clean(2:end);
RR_n_1 = RR_clean(1:end-1);

% SD1 and SD2
diff_RR = RR_n - RR_n_1;
sum_RR  = RR_n + RR_n_1;
SD1 = std(diff_RR)/sqrt(2);
SD2 = std(sum_RR)/sqrt(2);
SD1_SD2_ratio = SD1 / SD2;

fprintf('--- Poincaré Metrics (Full Data) ---\n');
fprintf('SD1 = %.4f s\n', SD1);
fprintf('SD2 = %.4f s\n', SD2);
fprintf('SD1/SD2 ratio = %.4f\n', SD1_SD2_ratio);

mean_x = mean(RR_n_1);
mean_y = mean(RR_n);

figure('Name','Poincare Full','Position',[100 100 700 700]);
hold on; grid on; axis equal;

% Scatter RR points
plot(RR_n_1, RR_n, 'o', 'MarkerFaceColor', [0.2 0.6 1], 'MarkerEdgeColor', 'k');

% Ellipse
theta = linspace(0,2*pi,200);
ellipse_x = SD2*cos(theta)/sqrt(2) - SD1*sin(theta)/sqrt(2) + mean_x;
ellipse_y = SD2*cos(theta)/sqrt(2) + SD1*sin(theta)/sqrt(2) + mean_y;
fill(ellipse_x, ellipse_y, [1 0.6 0.4], 'FaceAlpha',0.2, 'EdgeColor','r', 'LineWidth',1.5);

% Arrows SD1 and SD2
quiver(mean_x, mean_y, SD2, SD2, 'k','LineWidth',1.5,'MaxHeadSize',0.5);
quiver(mean_x, mean_y, -SD2, -SD2, 'k','LineWidth',1.5,'MaxHeadSize',0.5);
quiver(mean_x, mean_y, SD1, -SD1, 'g','LineWidth',1.5,'MaxHeadSize',0.5);
quiver(mean_x, mean_y, -SD1, SD1, 'g','LineWidth',1.5,'MaxHeadSize',0.5);

xlabel('RR[n-1] (s)');
ylabel('RR[n] (s)');
title('Poincaré Plot of RR Intervals (Full Data) with SD1/SD2');
legend('RR Points','SD1/SD2 Ellipse','SD2 Axis','SD1 Axis','Location','best');
hold off;

%% === Poincaré Zoomed (1-minute Interval) ===
zoom_start_time = 120; % 2 min
zoom_end_time   = 180; % 3 min
idx_zoom = R_locs/fs >= zoom_start_time & R_locs/fs <= zoom_end_time;

RR_zoom = diff(t(R_locs(idx_zoom)));
RR_zoom_clean = remove_outliers(RR_zoom);

RR_n   = RR_zoom_clean(2:end);
RR_n_1 = RR_zoom_clean(1:end-1);

% SD1 and SD2
diff_RR = RR_n - RR_n_1;
sum_RR  = RR_n + RR_n_1;
SD1 = std(diff_RR)/sqrt(2);
SD2 = std(sum_RR)/sqrt(2);
SD1_SD2_ratio = SD1 / SD2;

fprintf('--- Poincaré Metrics (Zoomed 2-3 min) ---\n');
fprintf('SD1 = %.4f s\n', SD1);
fprintf('SD2 = %.4f s\n', SD2);
fprintf('SD1/SD2 ratio = %.4f\n', SD1_SD2_ratio);

mean_x = mean(RR_n_1);
mean_y = mean(RR_n);

figure('Name','Poincare Zoomed','Position',[100 100 700 700]);
hold on; grid on; axis equal;

plot(RR_n_1, RR_n, 'o', 'MarkerFaceColor', [0.2 0.6 1], 'MarkerEdgeColor', 'k');

% Ellipse
ellipse_x = SD2*cos(theta)/sqrt(2) - SD1*sin(theta)/sqrt(2) + mean_x;
ellipse_y = SD2*cos(theta)/sqrt(2) + SD1*sin(theta)/sqrt(2) + mean_y;
fill(ellipse_x, ellipse_y, [1 0.6 0.4], 'FaceAlpha',0.2, 'EdgeColor','r', 'LineWidth',1.5);

% Arrows SD1 and SD2
quiver(mean_x, mean_y, SD2, SD2, 'k','LineWidth',1.5,'MaxHeadSize',0.5);
quiver(mean_x, mean_y, -SD2, -SD2, 'k','LineWidth',1.5,'MaxHeadSize',0.5);
quiver(mean_x, mean_y, SD1, -SD1, 'g','LineWidth',1.5,'MaxHeadSize',0.5);
quiver(mean_x, mean_y, -SD1, SD1, 'g','LineWidth',1.5,'MaxHeadSize',0.5);

xlabel('RR[n-1] (s)');
ylabel('RR[n] (s)');
title('Poincaré Plot of RR Intervals (Zoomed 2-3 min) with SD1/SD2');
legend('RR Points','SD1/SD2 Ellipse','SD2 Axis','SD1 Axis','Location','best');
hold off;
