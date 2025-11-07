%% ============================================================
%  ECG SIGNAL ANALYSIS - FILTERING, FREQUENCY RESPONSE & WAVE DETECTION
%  Application to Cardiovascular Fitness
%  Description:
%     Loads an ECG text file, filters it (0.5–40 Hz Butterworth),
%     plots the filter frequency responses,
%     detects P–QRS–T waves, and compares filtered leads I–II–III.
%% ============================================================

clc; clear; close all;

%% === Data file source ===
filename = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';

fid = fopen(filename, 'r');
if fid == -1
    error('File could not be opened. Check the path.');
end

% Skip header lines until "# EndOfHeader"
while true
    line = fgetl(fid);
    if ~ischar(line)
        error('"# EndOfHeader" not found in the file.');
    end
    if contains(line, '# EndOfHeader')
        break;
    end
end

% Read numerical data
data = textscan(fid, '%f%f%f%f%f%f', 'Delimiter', '\t', 'CollectOutput', true);
fclose(fid);
data = data{1};

fs = 100;  % Sampling frequency [Hz]

%% === Select ECG channel (A2 = Lead II) ===
ECG = data(:,6);
t = (0:length(ECG)-1)/fs;

%% === Plot raw ECG ===
figure;
plot(t, ECG, 'b');
xlabel('Time (s)'); ylabel('Amplitude (V)');
title('ECG - Lead II (Raw)');
grid on;

figure;
plot(t(1:5*fs), ECG(1:5*fs), 'r');
xlabel('Time (s)'); ylabel('Amplitude (V)');
title('Zoom - First 5 seconds of ECG');
grid on;

%% === Filtering ===
% High-pass filter (<0.5 Hz)
[b_hp, a_hp] = butter(2, 0.5/(fs/2), 'high');
ECG_hp = filtfilt(b_hp, a_hp, ECG);

% Low-pass filter (>40 Hz)
[b_lp, a_lp] = butter(2, 40/(fs/2), 'low');
ECG_filtered = filtfilt(b_lp, a_lp, ECG_hp);

% Convert to mV
ECG_mV = ECG * 1000;
ECG_filtered_mV = ECG_filtered * 1000;

%% === Frequency response of filters ===
figure;
freqz(b_hp, a_hp, 1024, fs);
title('High-pass Butterworth filter (0.5 Hz)');
grid on;

figure;
freqz(b_lp, a_lp, 1024, fs);
title('Low-pass Butterworth filter (40 Hz)');
grid on;

%% === Raw vs Filtered comparison ===
figure;
plot(t, ECG_mV, 'b'); hold on;
plot(t, ECG_filtered_mV, 'r');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('ECG - Raw (blue) vs Filtered (red)');
legend('Raw','Filtered');
grid on;

figure;
plot(t(1:5*fs), ECG_mV(1:5*fs), 'b'); hold on;
plot(t(1:5*fs), ECG_filtered_mV(1:5*fs), 'r');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('Zoom - First 5 seconds: Raw vs Filtered');
legend('Raw','Filtered');
grid on;

%% === QRS detection ===
minPeakHeight = 0.5 * max(ECG_filtered_mV);
minPeakDistance = 0.6 * fs;
[~, R_locs] = findpeaks(ECG_filtered_mV, ...
    'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDistance);

%% === P and T wave detection ===
P_locs = []; T_locs = [];
window_samples = round(0.2 * fs); % 200 ms

for i = 1:length(R_locs)
    if R_locs(i)-window_samples > 0
        [~, p_idx] = max(ECG_filtered_mV(R_locs(i)-window_samples:R_locs(i)-10));
        P_locs = [P_locs; R_locs(i)-window_samples-1 + p_idx];
    end
    if R_locs(i)+10+window_samples <= length(ECG_filtered_mV)
        [~, t_idx] = max(ECG_filtered_mV(R_locs(i)+10:R_locs(i)+window_samples));
        T_locs = [T_locs; R_locs(i)+10-1 + t_idx];
    end
end

%% === Plot detected waves (entire signal) ===
figure;
plot(t, ECG_filtered_mV, 'b'); hold on;
plot(t(R_locs), ECG_filtered_mV(R_locs), 'ro','MarkerFaceColor','r'); % R
plot(t(P_locs), ECG_filtered_mV(P_locs), 'go','MarkerFaceColor','g'); % P
plot(t(T_locs), ECG_filtered_mV(T_locs), 'mo','MarkerFaceColor','m'); % T
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('Filtered ECG with detected P, QRS, and T waves');
legend('ECG','R peaks','P waves','T waves');
grid on;

%% === Zoom (first 5 seconds) ===
zoom_window = 5; % seconds
idx = t <= zoom_window; % logical index

R_zoom = R_locs(R_locs <= find(idx,1,'last'));
P_zoom = P_locs(P_locs > 0 & P_locs <= find(idx,1,'last'));
T_zoom = T_locs(T_locs > 0 & T_locs <= find(idx,1,'last'));

figure;
plot(t(idx), ECG_filtered_mV(idx), 'b'); hold on;

if ~isempty(R_zoom)
    plot(t(R_zoom), ECG_filtered_mV(R_zoom), 'ro','MarkerFaceColor','r');
end
if ~isempty(P_zoom)
    plot(t(P_zoom), ECG_filtered_mV(P_zoom), 'go','MarkerFaceColor','g');
end
if ~isempty(T_zoom)
    plot(t(T_zoom), ECG_filtered_mV(T_zoom), 'mo','MarkerFaceColor','m');
end

xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('Zoom - First 5 seconds with detected P, QRS, T waves');
legend('ECG','R peaks','P waves','T waves');
grid on;

