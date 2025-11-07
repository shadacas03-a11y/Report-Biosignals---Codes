%% === ECG SIGNAL ANALYSIS (3 LEADS) ===
clear; close all; clc;

% === Data file sources ===
file1 = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';
file2 = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD2_resting_converted.txt';
file3 = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD3_resting_converted.txt';

% === Helper function to read ECG file ===
function ECG = readECGfile(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    % Skip header until "# EndOfHeader"
    while true
        line = fgetl(fid);
        if ~ischar(line)
            error('# EndOfHeader not found in %s', filename);
        end
        if contains(line, '# EndOfHeader')
            break;
        end
    end
    % Read numerical data
    data = textscan(fid, '%f%f%f%f%f%f', 'Delimiter', '\t', 'CollectOutput', true);
    fclose(fid);
    data = data{1};
    ECG = data(:,6); % Use channel 6 (Lead II format)
end

% === Load all 3 leads ===
LeadI  = readECGfile(file1);
LeadII = readECGfile(file2);
LeadIII = readECGfile(file3);

fs = 100; % Sampling frequency [Hz]
t = (0:length(LeadII)-1)/fs;

% === FILTER DESIGN ===
[b_hp, a_hp] = butter(2, 0.5/(fs/2), 'high');
[b_lp, a_lp] = butter(2, 40/(fs/2), 'low');

% === FREQUENCY RESPONSE OF FILTERS ===
figure;
freqz(b_hp, a_hp, 1024, fs);
title('High-pass Butterworth filter (0.5 Hz)');
grid on;

figure;
freqz(b_lp, a_lp, 1024, fs);
title('Low-pass Butterworth filter (40 Hz)');
grid on;

% === APPLY FILTERS ===
filt = @(x) filtfilt(b_lp, a_lp, filtfilt(b_hp, a_hp, x)) * 1000; % filtered in mV
ECG_I_filt   = filt(LeadI);
ECG_II_filt  = filt(LeadII);
ECG_III_filt = filt(LeadIII);

% === PLOT FILTERED LEADS ===
figure('Position',[100 100 1200 400]);
plot(t, ECG_I_filt, 'r'); hold on;
plot(t, ECG_II_filt, 'b');
plot(t, ECG_III_filt, 'g');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('Filtered ECG signals from Leads I, II, and III');
legend('Lead I', 'Lead II', 'Lead III');
grid on;

% Zoom on first 5 seconds
figure('Position',[100 100 1200 400]);
idx_zoom = t <= 5;
plot(t(idx_zoom), ECG_I_filt(idx_zoom), 'r'); hold on;
plot(t(idx_zoom), ECG_II_filt(idx_zoom), 'b');
plot(t(idx_zoom), ECG_III_filt(idx_zoom), 'g');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('Zoom - Filtered ECG (first 5 seconds)');
legend('Lead I', 'Lead II', 'Lead III');
grid on;

% === QRS DETECTION (Lead II only) ===
ECG_filtered_mV = ECG_II_filt;
minPeakHeight = 0.5 * max(ECG_filtered_mV);
minPeakDistance = 0.6 * fs;
[~, R_locs] = findpeaks(ECG_filtered_mV, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDistance);

% === P and T waves (approximation) ===
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

% === PLOT DETECTED WAVES (Lead II) ===
figure;
plot(t, ECG_filtered_mV, 'b'); hold on;
plot(t(R_locs), ECG_filtered_mV(R_locs), 'ro','MarkerFaceColor','r'); % R
plot(t(P_locs), ECG_filtered_mV(P_locs), 'go','MarkerFaceColor','g'); % P
plot(t(T_locs), ECG_filtered_mV(T_locs), 'mo','MarkerFaceColor','m'); % T
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('Filtered ECG (Lead II) with P, QRS, T waves');
legend('ECG','R peaks','P waves','T waves');
grid on;

% === ZOOM on first 5 seconds (Lead II) ===
zoom_window = 5; % seconds
idx = t <= zoom_window;
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
title('Zoom - First 5 seconds of filtered ECG (Lead II)');
legend('ECG','R peaks','P waves','T waves');
grid on;
