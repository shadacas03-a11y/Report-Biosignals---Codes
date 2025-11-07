%% === Comparative ECG Analysis: Student A vs Student B ===
clc; clear; close all;

%% === PARAMETERS ===
fs = 100; % ECG sampling frequency [Hz]
fs_resample = 4; % resample frequency for RR intervals
VLF_band = [0 0.04];
LF_band  = [0.04 0.15];
HF_band  = [0.15 0.4];

%% === FUNCTION TO REMOVE OUTLIERS ===
remove_outliers = @(x) x(x >= prctile(x,25)-1.5*(prctile(x,75)-prctile(x,25)) & ...
                          x <= prctile(x,75)+1.5*(prctile(x,75)-prctile(x,25)));

%% === Load Student A ECG ===
filename_A = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';
data_A = readmatrix(filename_A);
ECG_A = data_A(:,6);
t_A = (0:length(ECG_A)-1)/fs;

% Filter
[b_hp,a_hp] = butter(2,0.5/(fs/2),'high');
[b_lp,a_lp] = butter(4,35/(fs/2),'low');
ECG_A_filt = filtfilt(b_lp,a_lp,filtfilt(b_hp,a_hp,ECG_A))*1000;
ECG_A_smooth = smoothdata(ECG_A_filt,'gaussian',round(0.015*fs));

% R-Peak Detection
minPeakHeight = 0.25*max(ECG_A_smooth);
minPeakDistance = 0.5*fs;
[R_peaks_A, R_locs_A] = findpeaks(ECG_A_smooth,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);

% RR intervals
RR_A = diff(t_A(R_locs_A));
RR_clean_A = remove_outliers(RR_A);

% Time-domain metrics
mean_RR_A = mean(RR_clean_A);
SDNN_A = std(RR_clean_A);
RMSSD_A = sqrt(mean(diff(RR_clean_A).^2))*1000;
avg_HR_A = mean(60 ./ RR_clean_A);

% Spectral Analysis
t_uniform_A = 0:1/fs_resample:sum(RR_clean_A);
RR_resampled_A = interp1(cumsum(RR_clean_A), RR_clean_A, t_uniform_A,'pchip');
RR_detrended_A = detrend(RR_resampled_A);
[Pxx_A, f_A] = pwelch(RR_detrended_A,[],[],[],fs_resample);

VLF_power_A = bandpower(Pxx_A,f_A,VLF_band,'psd');
LF_power_A  = bandpower(Pxx_A,f_A,LF_band,'psd');
HF_power_A  = bandpower(Pxx_A,f_A,HF_band,'psd');
LF_HF_ratio_A = LF_power_A/HF_power_A;

% Poincaré metrics
RR_n_A = RR_clean_A(2:end);
RR_n1_A = RR_clean_A(1:end-1);
SD1_A = std(RR_n_A - RR_n1_A)/sqrt(2);
SD2_A = std(RR_n_A + RR_n1_A)/sqrt(2);
SD1_SD2_ratio_A = SD1_A/SD2_A;

%% === Load Student B ECG ===
filename_B = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Rodrigo_LEAD1_resting_converted.txt';
data_B = readmatrix(filename_B);
ECG_B = data_B(:,6);
t_B = (0:length(ECG_B)-1)/fs;

% Filter
ECG_B_filt = filtfilt(b_lp,a_lp,filtfilt(b_hp,a_hp,ECG_B))*1000;
ECG_B_smooth = smoothdata(ECG_B_filt,'gaussian',round(0.015*fs));

% R-Peak Detection
[R_peaks_B, R_locs_B] = findpeaks(ECG_B_smooth,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);

% RR intervals
RR_B = diff(t_B(R_locs_B));
RR_clean_B = remove_outliers(RR_B);

% Time-domain metrics
mean_RR_B = mean(RR_clean_B);
SDNN_B = std(RR_clean_B);
RMSSD_B = sqrt(mean(diff(RR_clean_B).^2))*1000;
avg_HR_B = mean(60 ./ RR_clean_B);

% Spectral Analysis
t_uniform_B = 0:1/fs_resample:sum(RR_clean_B);
RR_resampled_B = interp1(cumsum(RR_clean_B), RR_clean_B, t_uniform_B,'pchip');
RR_detrended_B = detrend(RR_resampled_B);
[Pxx_B, f_B] = pwelch(RR_detrended_B,[],[],[],fs_resample);

VLF_power_B = bandpower(Pxx_B,f_B,VLF_band,'psd');
LF_power_B  = bandpower(Pxx_B,f_B,LF_band,'psd');
HF_power_B  = bandpower(Pxx_B,f_B,HF_band,'psd');
LF_HF_ratio_B = LF_power_B/HF_power_B;

% Poincaré metrics
RR_n_B = RR_clean_B(2:end);
RR_n1_B = RR_clean_B(1:end-1);
SD1_B = std(RR_n_B - RR_n1_B)/sqrt(2);
SD2_B = std(RR_n_B + RR_n1_B)/sqrt(2);
SD1_SD2_ratio_B = SD1_B/SD2_B;

%% === Comparative Table ===
Metric = {'Mean RR (s)'; 'SDNN (s)'; 'RMSSD (ms)'; 'Avg HR (bpm)'; ...
          'VLF Power (s^2)'; 'LF Power (s^2)'; 'HF Power (s^2)'; 'LF/HF'; ...
          'SD1 (s)'; 'SD2 (s)'; 'SD1/SD2'};
Student_A = [mean_RR_A; SDNN_A; RMSSD_A; avg_HR_A; VLF_power_A; LF_power_A; HF_power_A; LF_HF_ratio_A; SD1_A; SD2_A; SD1_SD2_ratio_A];
Student_B = [mean_RR_B; SDNN_B; RMSSD_B; avg_HR_B; VLF_power_B; LF_power_B; HF_power_B; LF_HF_ratio_B; SD1_B; SD2_B; SD1_SD2_ratio_B];

T = table(Metric, Student_A, Student_B)
disp('--- Comparative Table ---');

%% === Plot Comparison Figures ===

% 1. RR Evolution
figure('Name','RR Interval Evolution Comparison','Position',[100 100 1200 400]);
plot(cumsum(RR_clean_A), 60./RR_clean_A,'-o','LineWidth',1.2); hold on;
plot(cumsum(RR_clean_B), 60./RR_clean_B,'-s','LineWidth',1.2);
xlabel('Time (s)'); ylabel('Heart Rate (bpm)');
title('Heart Rate Evolution: Student A vs Student B');
legend('Student A','Student B'); grid on;

% 2. RR Resampled
figure('Name','RR Resampled Comparison','Position',[100 100 1200 400]);
plot(t_uniform_A, RR_resampled_A,'-o','LineWidth',1.2); hold on;
plot(t_uniform_B, RR_resampled_B,'-s','LineWidth',1.2);
xlabel('Time (s)'); ylabel('RR Interval (s)');
title('RR Intervals Resampled at 4 Hz: Student A vs Student B');
legend('Student A','Student B'); grid on;

% 3. PSD Comparison
figure('Name','PSD Comparison','Position',[100 100 1200 400]);
hold on; grid on;

y_max = max([Pxx_A; Pxx_B])*1.1;

% Colored bands
patch([VLF_band(1) VLF_band(2) VLF_band(2) VLF_band(1)], [0 0 y_max y_max],[0.9 0.9 0.9],'FaceAlpha',0.3,'EdgeColor','none');
patch([LF_band(1) LF_band(2) LF_band(2) LF_band(1)], [0 0 y_max y_max],[0.6 0.8 1],'FaceAlpha',0.3,'EdgeColor','none');
patch([HF_band(1) HF_band(2) HF_band(2) HF_band(1)], [0 0 y_max y_max],[0.8 1 0.6],'FaceAlpha',0.3,'EdgeColor','none');

plot(f_A, Pxx_A,'k','LineWidth',1.5);
plot(f_B, Pxx_B,'r','LineWidth',1.5);

xlabel('Frequency (Hz)'); ylabel('Power Spectral Density (s^2/Hz)');
title('PSD Comparison of RR Intervals: Student A vs Student B');
legend('VLF','LF','HF','Student A','Student B'); xlim([0 0.5]); ylim([0 y_max]);

% 4. Poincaré plots side-by-side
figure('Name','Poincaré Plot Comparison','Position',[100 100 1200 500]);
theta = linspace(0,2*pi,200);

% Student A
subplot(1,2,1); hold on; grid on; axis equal;
plot(RR_n1_A, RR_n_A,'o','MarkerFaceColor',[0.2 0.6 1],'MarkerEdgeColor','k');
ellipse_x = SD2_A*cos(theta)/sqrt(2) - SD1_A*sin(theta)/sqrt(2) + mean(RR_n1_A);
ellipse_y = SD2_A*cos(theta)/sqrt(2) + SD1_A*sin(theta)/sqrt(2) + mean(RR_n_A);
fill(ellipse_x, ellipse_y,[1 0.6 0.4],'FaceAlpha',0.2,'EdgeColor','r','LineWidth',1.5);
quiver(mean(RR_n1_A), mean(RR_n_A), SD2_A, SD2_A,'k','LineWidth',1.5,'MaxHeadSize',0.5);
quiver(mean(RR_n1_A), mean(RR_n_A), SD1_A, -SD1_A,'g','LineWidth',1.5,'MaxHeadSize',0.5);
xlabel('RR[n-1] (s)'); ylabel('RR[n] (s)'); title('Student A'); legend('RR Points','SD1/SD2 Ellipse','SD2','SD1');

% Student B
subplot(1,2,2); hold on; grid on; axis equal;
plot(RR_n1_B, RR_n_B,'o','MarkerFaceColor',[0.2 0.6 1],'MarkerEdgeColor','k');
ellipse_x_B = SD2_B*cos(theta)/sqrt(2) - SD1_B*sin(theta)/sqrt(2) + mean(RR_n1_B);
ellipse_y_B = SD2_B*cos(theta)/sqrt(2) + SD1_B*sin(theta)/sqrt(2) + mean(RR_n_B);
fill(ellipse_x_B, ellipse_y_B,[1 0.6 0.4],'FaceAlpha',0.2,'EdgeColor','r','LineWidth',1.5);
quiver(mean(RR_n1_B), mean(RR_n_B), SD2_B, SD2_B,'k','LineWidth',1.5,'MaxHeadSize',0.5);
quiver(mean(RR_n1_B), mean(RR_n_B), SD1_B, -SD1_B,'g','LineWidth',1.5,'MaxHeadSize',0.5);
xlabel('RR[n-1] (s)'); ylabel('RR[n] (s)'); title('Student B'); legend('RR Points','SD1/SD2 Ellipse','SD2','SD1');
