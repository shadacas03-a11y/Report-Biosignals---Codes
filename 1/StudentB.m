%% === Student B: ECG Analysis with Full PSD Colors and Legend ===
clc; clear; close all;

%% Load ECG
filename_B = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Rodrigo_LEAD1_resting_converted.txt';
data_B = readmatrix(filename_B);
fs = 100;                   % Sampling frequency [Hz]
ECG_B = data_B(:,6);        % Lead III
t_B = (0:length(ECG_B)-1)/fs;

%% Filter ECG
[b_hp,a_hp] = butter(2,0.5/(fs/2),'high');
[b_lp,a_lp] = butter(4,35/(fs/2),'low');
ECG_B_filt = filtfilt(b_lp,a_lp,filtfilt(b_hp,a_hp,ECG_B))*1000;
ECG_B_smooth = smoothdata(ECG_B_filt,'gaussian',round(0.015*fs));

%% Detect R-peaks
minPeakHeight = 0.25*max(ECG_B_smooth);
minPeakDistance = 0.5*fs;
[R_peaks_B, R_locs_B] = findpeaks(ECG_B_smooth,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);

%% RR Intervals & Remove Outliers
RR_B = diff(t_B(R_locs_B));
remove_outliers = @(x) x(x >= prctile(x,25)-1.5*(prctile(x,75)-prctile(x,25)) & ...
                          x <= prctile(x,75)+1.5*(prctile(x,75)-prctile(x,25)));
RR_clean_B = remove_outliers(RR_B);

%% Time-Domain Metrics
mean_RR_B = mean(RR_clean_B);
SDNN_B = std(RR_clean_B);
RMSSD_B = sqrt(mean(diff(RR_clean_B).^2))*1000; % ms
avg_HR_B = mean(60 ./ RR_clean_B);

fprintf('--- Student B Time-Domain Metrics ---\n');
fprintf('Mean RR = %.3f s\nSDNN = %.4f s\nRMSSD = %.2f ms\nAverage HR = %.1f bpm\n', ...
        mean_RR_B, SDNN_B, RMSSD_B, avg_HR_B);

%% Spectral Analysis
fs_resample = 4; % resampling at 4 Hz
t_uniform_B = 0:1/fs_resample:sum(RR_clean_B);
RR_resampled_B = interp1(cumsum(RR_clean_B), RR_clean_B, t_uniform_B, 'pchip');
RR_detrended_B = detrend(RR_resampled_B);

[Pxx_B,f_B] = pwelch(RR_detrended_B,[],[],[],fs_resample);

% Define frequency bands
VLF_band = [0 0.04];
LF_band  = [0.04 0.15];
HF_band  = [0.15 0.4];

VLF_power_B = bandpower(Pxx_B,f_B,VLF_band,'psd');
LF_power_B  = bandpower(Pxx_B,f_B,LF_band,'psd');
HF_power_B  = bandpower(Pxx_B,f_B,HF_band,'psd');
LF_HF_ratio_B = LF_power_B/HF_power_B;

fprintf('--- Student B Spectral Metrics ---\n');
fprintf('VLF Power = %.4f s^2\nLF Power = %.4f s^2\nHF Power = %.4f s^2\nLF/HF Ratio = %.2f\n', ...
        VLF_power_B, LF_power_B, HF_power_B, LF_HF_ratio_B);

%% Poincaré Plot Metrics
RR_n_B  = RR_clean_B(2:end);
RR_n1_B = RR_clean_B(1:end-1);

SD1_B = std(RR_n_B - RR_n1_B)/sqrt(2);
SD2_B = std(RR_n_B + RR_n1_B)/sqrt(2);
SD1_SD2_ratio_B = SD1_B/SD2_B;

fprintf('--- Student B Poincaré Metrics ---\n');
fprintf('SD1 = %.4f s\nSD2 = %.4f s\nSD1/SD2 ratio = %.4f\n', SD1_B, SD2_B, SD1_SD2_ratio_B);

%% === RR Interval Evolution Plot ===
figure('Name','Student B - RR Intervals','Position',[100 100 800 400]);
t_RR_B = cumsum(RR_clean_B);
HR_series_B = 60 ./ RR_clean_B;
plot(t_RR_B, HR_series_B, '-o','LineWidth',1.2,'MarkerSize',4);
xlabel('Time (s)'); ylabel('Heart Rate (bpm)');
title('Student B - Heart Rate Evolution Over Time');
grid on;

%% === RR Resampled Plot ===
figure('Name','Student B - RR Resampled 4Hz','Position',[100 100 800 400]);
plot(t_uniform_B, RR_resampled_B, '-o','LineWidth',1.2,'MarkerSize',4);
xlabel('Time (s)'); ylabel('RR Interval (s)');
title('Student B - RR Intervals Resampled at 4 Hz');
grid on;

%% === PSD with Colored Bands and Legend ===
figure('Name','Student B - PSD','Position',[100 100 800 400]);
hold on; grid on;

y_max = max(Pxx_B)*1.1; % margin for patches

% Draw colored bands
hVLF = patch([VLF_band(1) VLF_band(2) VLF_band(2) VLF_band(1)], [0 0 y_max y_max], [0.9 0.9 0.9], 'FaceAlpha',0.3, 'EdgeColor','none');
hLF  = patch([LF_band(1) LF_band(2) LF_band(2) LF_band(1)], [0 0 y_max y_max], [0.6 0.8 1], 'FaceAlpha',0.3, 'EdgeColor','none');
hHF  = patch([HF_band(1) HF_band(2) HF_band(2) HF_band(1)], [0 0 y_max y_max], [0.8 1 0.6], 'FaceAlpha',0.3, 'EdgeColor','none');

% PSD line
hPSD = plot(f_B,Pxx_B,'k','LineWidth',1.5);

% Add text labels
text(mean(VLF_band), 0.95*y_max, 'VLF','HorizontalAlignment','center','FontWeight','bold');
text(mean(LF_band), 0.95*y_max, 'LF','HorizontalAlignment','center','FontWeight','bold');
text(mean(HF_band), 0.95*y_max, 'HF','HorizontalAlignment','center','FontWeight','bold');

xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (s^2/Hz)');
title('Student B - PSD of RR Intervals with VLF/LF/HF Bands');
xlim([0 0.5]);
ylim([0 y_max]);

% Legend
legend([hVLF, hLF, hHF, hPSD], {'VLF band','LF band','HF band','PSD'}, 'Location','northeast');

hold off;

%% === Poincaré Plot ===
mean_x_B = mean(RR_n1_B);
mean_y_B = mean(RR_n_B);

figure('Name','Student B - Poincare Plot','Position',[100 100 700 700]);
hold on; grid on; axis equal;

% Scatter plot
plot(RR_n1_B, RR_n_B, 'o', 'MarkerFaceColor',[0.2 0.6 1], 'MarkerEdgeColor','k');

% SD1/SD2 ellipse
theta = linspace(0,2*pi,200);
ellipse_x_B = SD2_B*cos(theta)/sqrt(2) - SD1_B*sin(theta)/sqrt(2) + mean_x_B;
ellipse_y_B = SD2_B*cos(theta)/sqrt(2) + SD1_B*sin(theta)/sqrt(2) + mean_y_B;
fill(ellipse_x_B, ellipse_y_B, [1 0.6 0.4], 'FaceAlpha',0.2,'EdgeColor','r','LineWidth',1.5);

% Axes arrows
quiver(mean_x_B, mean_y_B, SD2_B, SD2_B, 'k','LineWidth',1.5,'MaxHeadSize',0.5); % SD2
quiver(mean_x_B, mean_y_B, -SD2_B, -SD2_B, 'k','LineWidth',1.5,'MaxHeadSize',0.5);
quiver(mean_x_B, mean_y_B, SD1_B, -SD1_B, 'g','LineWidth',1.5,'MaxHeadSize',0.5); % SD1
quiver(mean_x_B, mean_y_B, -SD1_B, SD1_B, 'g','LineWidth',1.5,'MaxHeadSize',0.5);

xlabel('RR[n-1] (s)'); ylabel('RR[n] (s)');
title('Student B - Poincaré Plot with SD1/SD2 Geometry');
legend('RR Points','SD1/SD2 Ellipse','SD2 Axis','SD1 Axis','Location','best');

hold off;
