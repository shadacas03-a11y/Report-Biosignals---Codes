%% === Student A & B ECG Analysis with RR Interval Comparison ===
clc; clear; close all;

%% --- Load Student A ECG ---
filename_A = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';
data_A = readmatrix(filename_A);
fs = 100; % Sampling frequency [Hz]
ECG_A = data_A(:,6);
t_A = (0:length(ECG_A)-1)/fs;

%% --- Filter ECG A ---
[b_hp,a_hp] = butter(2,0.5/(fs/2),'high');
[b_lp,a_lp] = butter(4,35/(fs/2),'low');
ECG_A_filt = filtfilt(b_lp,a_lp,filtfilt(b_hp,a_hp,ECG_A))*1000;
ECG_A_smooth = smoothdata(ECG_A_filt,'gaussian',round(0.015*fs));

%% --- R-Peak Detection A ---
minPeakHeight = 0.25*max(ECG_A_smooth);
minPeakDistance = 0.5*fs;
[R_peaks_A,R_locs_A] = findpeaks(ECG_A_smooth,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);

%% --- RR Intervals A & Remove Outliers ---
RR_A = diff(t_A(R_locs_A));
remove_outliers = @(x) x(x >= prctile(x,25)-1.5*(prctile(x,75)-prctile(x,25)) & ...
                          x <= prctile(x,75)+1.5*(prctile(x,75)-prctile(x,25)));
RR_clean_A = remove_outliers(RR_A);

%% --- Time-Domain Metrics A ---
mean_RR_A = mean(RR_clean_A);
SDNN_A = std(RR_clean_A);
RMSSD_A = sqrt(mean(diff(RR_clean_A).^2))*1000; % ms
avg_HR_A = mean(60 ./ RR_clean_A);

fprintf('--- Student A Time-Domain Metrics ---\n');
fprintf('Mean RR = %.3f s\nSDNN = %.4f s\nRMSSD = %.2f ms\nAverage HR = %.1f bpm\n', ...
        mean_RR_A, SDNN_A, RMSSD_A, avg_HR_A);

%% --- Load Student B ECG ---
filename_B = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Rodrigo_LEAD1_resting_converted.txt';
data_B = readmatrix(filename_B);
ECG_B = data_B(:,6);
t_B = (0:length(ECG_B)-1)/fs;

%% --- Filter ECG B ---
ECG_B_filt = filtfilt(b_lp,a_lp,filtfilt(b_hp,a_hp,ECG_B))*1000;
ECG_B_smooth = smoothdata(ECG_B_filt,'gaussian',round(0.015*fs));

%% --- R-Peak Detection B ---
[R_peaks_B,R_locs_B] = findpeaks(ECG_B_smooth,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);

%% --- RR Intervals B & Remove Outliers ---
RR_B = diff(t_B(R_locs_B));
RR_clean_B = remove_outliers(RR_B);

%% --- Time-Domain Metrics B ---
mean_RR_B = mean(RR_clean_B);
SDNN_B = std(RR_clean_B);
RMSSD_B = sqrt(mean(diff(RR_clean_B).^2))*1000; % ms
avg_HR_B = mean(60 ./ RR_clean_B);

fprintf('--- Student B Time-Domain Metrics ---\n');
fprintf('Mean RR = %.3f s\nSDNN = %.4f s\nRMSSD = %.2f ms\nAverage HR = %.1f bpm\n', ...
        mean_RR_B, SDNN_B, RMSSD_B, avg_HR_B);

%% --- PSD Analysis for Both Students ---
fs_resample = 4;

% Student A
t_uniform_A = 0:1/fs_resample:sum(RR_clean_A);
RR_resampled_A = interp1(cumsum(RR_clean_A), RR_clean_A, t_uniform_A, 'pchip');
RR_detrended_A = detrend(RR_resampled_A);
[Pxx_A,f_A] = pwelch(RR_detrended_A,[],[],[],fs_resample);

% Student B
t_uniform_B = 0:1/fs_resample:sum(RR_clean_B);
RR_resampled_B = interp1(cumsum(RR_clean_B), RR_clean_B, t_uniform_B, 'pchip');
RR_detrended_B = detrend(RR_resampled_B);
[Pxx_B,f_B] = pwelch(RR_detrended_B,[],[],[],fs_resample);

% Frequency bands
VLF_band = [0 0.04];
LF_band = [0.04 0.15];
HF_band = [0.15 0.4];

%% --- Compute Band Powers ---
VLF_A = bandpower(Pxx_A,f_A,VLF_band,'psd');
LF_A = bandpower(Pxx_A,f_A,LF_band,'psd');
HF_A = bandpower(Pxx_A,f_A,HF_band,'psd');
LF_HF_ratio_A = LF_A/HF_A;

VLF_B = bandpower(Pxx_B,f_B,VLF_band,'psd');
LF_B = bandpower(Pxx_B,f_B,LF_band,'psd');
HF_B = bandpower(Pxx_B,f_B,HF_band,'psd');
LF_HF_ratio_B = LF_B/HF_B;

fprintf('--- Student A Spectral Metrics ---\nVLF = %.4f, LF = %.4f, HF = %.4f, LF/HF = %.2f\n',VLF_A,LF_A,HF_A,LF_HF_ratio_A);
fprintf('--- Student B Spectral Metrics ---\nVLF = %.4f, LF = %.4f, HF = %.4f, LF/HF = %.2f\n',VLF_B,LF_B,HF_B,LF_HF_ratio_B);

%% --- Poincaré Plots ---
% Student A
RR_n_A = RR_clean_A(2:end); RR_n1_A = RR_clean_A(1:end-1);
SD1_A = std(RR_n_A - RR_n1_A)/sqrt(2);
SD2_A = std(RR_n_A + RR_n1_A)/sqrt(2);
SD1_SD2_ratio_A = SD1_A/SD2_A;

% Student B
RR_n_B = RR_clean_B(2:end); RR_n1_B = RR_clean_B(1:end-1);
SD1_B = std(RR_n_B - RR_n1_B)/sqrt(2);
SD2_B = std(RR_n_B + RR_n1_B)/sqrt(2);
SD1_SD2_ratio_B = SD1_B/SD2_B;

%% --- RR Interval Comparison (Statistical Test) ---
% Normality test
[hA,pA] = lillietest(RR_clean_A);
[hB,pB] = lillietest(RR_clean_B);

fprintf('--- Normality Test ---\n');
fprintf('Student A: h=%d, p=%.3f -> %s\n', hA, pA, ternary(hA==0,'Normal','Not Normal'));
fprintf('Student B: h=%d, p=%.3f -> %s\n', hB, pB, ternary(hB==0,'Normal','Not Normal'));

% Choose test
if hA==0 && hB==0
    [h_test,p_test] = ttest2(RR_clean_A, RR_clean_B);
    test_name = 'Two-sample t-test';
else
    [p_test,h_test] = ranksum(RR_clean_A, RR_clean_B);
    test_name = 'Mann-Whitney U test';
end

fprintf('--- Comparison of RR Intervals ---\n');
fprintf('Test used: %s\n', test_name);
fprintf('p-value = %.4f\n', p_test);
if p_test < 0.05
    fprintf('Conclusion: RR intervals differ significantly between Students A and B.\n');
else
    fprintf('Conclusion: No significant difference in RR intervals between Students A and B.\n');
end

%% --- Helper function ---
function out = ternary(cond, valTrue, valFalse)
    if cond, out = valTrue; else, out = valFalse; end
end
