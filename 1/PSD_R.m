%% === ECG P-Q-R-S-T Detection, Interval Analysis, and HRV Spectral Analysis ===
% -------------------------------------------------------
clc; clear; close all;

%% === CONFIGURATION ===
filename = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';
fs = 100; % ECG sampling frequency [Hz]
zoom_start = 176; zoom_end = 182; % segment for detection (adjust as needed)

%% === LOAD DATA ===
fid = fopen(filename,'r');
if fid == -1
    error('Cannot open file. Check the path.');
end
% Skip header until '# EndOfHeader'
while true
    line = fgetl(fid);
    if ~ischar(line), error('# EndOfHeader not found.'); end
    if contains(line,'# EndOfHeader'), break; end
end
data = textscan(fid,'%f%f%f%f%f%f','Delimiter','\t','CollectOutput',true);
fclose(fid);
data = data{1};
ECG = data(:,6);           
t = (0:length(ECG)-1)/fs;

%% === FILTERING ===
[b_hp,a_hp] = butter(2,0.5/(fs/2),'high');
[b_lp,a_lp] = butter(4,35/(fs/2),'low');
ECG_filt = filtfilt(b_lp,a_lp,filtfilt(b_hp,a_hp,ECG))*1000;  % convert to mV
ECG_smooth = smoothdata(ECG_filt,'gaussian',round(0.015*fs));

%% === SEGMENT SELECTION (ZOOM) ===
idx_zoom = t >= zoom_start & t <= zoom_end;
t_zoom = t(idx_zoom);
ECG_zoom = ECG_smooth(idx_zoom);
if isempty(ECG_zoom)
    error('Zoom interval is empty. Adjust zoom_start/zoom_end.');
end

%% === P-Q-R-S-T DETECTION ===
minPeakHeight = 0.25 * max(ECG_zoom);
minPeakDistance = 0.5 * fs; 
[R_peaks, R_locs] = findpeaks(ECG_zoom,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);
[~, min_locs] = findpeaks(-ECG_zoom);

P_locs = []; Q_locs = []; S_locs = []; T_locs = [];

for i = 1:length(R_locs)
    r = R_locs(i);
    % Q before R
    q_candidates = min_locs(min_locs < r);
    if ~isempty(q_candidates)
        q = q_candidates(end); Q_locs = [Q_locs q];
    else, q = NaN;
    end
    % S after R
    s_candidates = min_locs(min_locs > r);
    if ~isempty(s_candidates)
        s = s_candidates(1); S_locs = [S_locs s];
    else, s = NaN;
    end
    % P before Q
    if ~isnan(q)
        sr = max(q - round(0.25 * fs), 1):q;
        [pksP, locsP] = findpeaks(ECG_zoom(sr));
        if ~isempty(pksP)
            P_locs = [P_locs sr(locsP(end))];
        end
    end
    % T after S
    if ~isnan(s)
        sr = s:min(s + round(0.5 * fs), length(ECG_zoom));
        [pksT, locsT] = findpeaks(ECG_zoom(sr));
        if ~isempty(pksT)
            T_locs = [T_locs sr(locsT(1))];
        end
    end
end

%% === INTERVALS (seconds) ===
R_times = t_zoom(R_locs(~isnan(R_locs)));
Q_times = t_zoom(Q_locs(~isnan(Q_locs)));
S_times = t_zoom(S_locs(~isnan(S_locs)));
T_times = t_zoom(T_locs(~isnan(T_locs)));

RR_intervals = diff(R_times);
QQ_intervals = diff(Q_times);
SS_intervals = diff(S_times);
TT_intervals = diff(T_times);

%% === REMOVE OUTLIERS ===
remove_outliers = @(x) x(x >= prctile(x,25)-1.5*(prctile(x,75)-prctile(x,25)) & ...
                           x <= prctile(x,75)+1.5*(prctile(x,75)-prctile(x,25)));
RR_clean = remove_outliers(RR_intervals);
QQ_clean = remove_outliers(QQ_intervals);
SS_clean = remove_outliers(SS_intervals);
TT_clean = remove_outliers(TT_intervals);

%% === STATISTICS ===
intervals_data = {RR_clean, QQ_clean, SS_clean, TT_clean};
names = {'RR','QQ','SS','TT'};
fprintf('\n--- Interval Statistics ---\n');
for i = 1:length(intervals_data)
    d = intervals_data{i};
    if isempty(d)
        fprintf('%s-interval: No valid data available.\n', names{i});
        continue;
    end
    fprintf('%s-interval: mean = %.3f s, std = %.3f s, median = %.3f s\n', names{i}, mean(d), std(d), median(d));
end

%% === RESAMPLE RR AT 4 Hz ===
if ~isempty(RR_clean)
    t_RR = cumsum(RR_clean);
    fs_resample = 4; 
    t_uniform = t_RR(1):1/fs_resample:t_RR(end);
    RR_resampled = interp1(t_RR, RR_clean, t_uniform, 'pchip');

    figure('Position',[100 100 800 360]);
    plot(t_uniform, RR_resampled, '-o','LineWidth',1.2,'MarkerSize',4);
    xlabel('Time (s)'); ylabel('RR Interval (s)');
    title('RR Intervals Resampled at 4 Hz');
    grid on;
else
    error('RR_clean is empty: not enough valid RR intervals to proceed.');
end

%% === SPECTRAL ANALYSIS (Welch) ===
if exist('RR_resampled','var') && ~isempty(RR_resampled)
    RR_ms = RR_resampled * 1000; % to milliseconds
    RR_detrended = detrend(RR_ms);
    L = length(RR_detrended);
    if L < 16
        error('RR signal too short for spectral analysis.');
    end

    win_length = min(256, max(8, floor(L/2)));
    window = hamming(win_length);
    noverlap = round(0.5 * win_length);
    nfft = max(512, 2^nextpow2(L));
    fs_r = fs_resample;
    [Pxx,f] = pwelch(RR_detrended, window, noverlap, nfft, fs_r);

    % Frequency bands
    VLF_band = [0 0.04];
    LF_band  = [0.04 0.15];
    HF_band  = [0.15 0.4];

    % Power
    VLF_power = bandpower(Pxx,f,VLF_band,'psd');
    LF_power  = bandpower(Pxx,f,LF_band,'psd');
    HF_power  = bandpower(Pxx,f,HF_band,'psd');
    Total_power = bandpower(Pxx,f,[0 0.4],'psd');
    LF_HF_ratio = LF_power / HF_power;

    fprintf('\n=== Spectral Analysis ===\n');
    fprintf('VLF Power (0–0.04 Hz): %.4f ms^2\n', VLF_power);
    fprintf('LF  Power (0.04–0.15 Hz): %.4f ms^2\n', LF_power);
    fprintf('HF  Power (0.15–0.4 Hz): %.4f ms^2\n', HF_power);
    fprintf('Total Power (0–0.4 Hz): %.4f ms^2\n', Total_power);
    fprintf('LF/HF Ratio: %.2f\n', LF_HF_ratio);

    %% === PSD PLOT WITH COLORED BANDS AND CORRECT LEGEND ===
PSD_RR = figure('Position',[120 120 920 420]); hold on;

f = f(:); Pxx = Pxx(:);

% Define band indices
idxV = f >= VLF_band(1) & f <= VLF_band(2);
idxL = f >= LF_band(1) & f <= LF_band(2);
idxH = f >= HF_band(1) & f <= HF_band(2);

% Function to fill bands
plot_band = @(idx, color) fill([f(idx); flipud(f(idx))], [Pxx(idx); zeros(sum(idx),1)], color, ...
                               'FaceAlpha',0.7, 'EdgeColor','none');

% Draw bands and save handles for legend
hVLF = []; hLF = []; hHF = [];
if any(idxV), hVLF = plot_band(idxV, [1 0.6 0.8]); end   % pink
if any(idxL), hLF  = plot_band(idxL, [0.4 0.75 1]); end  % blue
if any(idxH), hHF  = plot_band(idxH, [1 0.8 0.5]); end   % orange

% Plot PSD curve and save handle
hPSD = plot(f, Pxx, 'k', 'LineWidth', 1.5);

% Correct legend using handles
legend_handles = [hPSD hVLF hLF hHF];
legend_labels = {'PSD','VLF','LF','HF'};
legend(legend_handles, legend_labels, 'Location','northeastoutside');

% Labels
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Power Spectral Density (ms^2/Hz)','FontSize',12);
title('PSD of RR Intervals (Welch Periodogram)','FontSize',14);
grid on;

% Optional: text labels above bands
ymax = max(Pxx); ypos = ymax * 0.78;
if any(idxV), text(mean(VLF_band), ypos, 'VLF', 'HorizontalAlignment','center','FontSize',12,'FontWeight','bold','Color',[0.6 0 0.4]); end
if any(idxL), text(mean(LF_band), ypos, 'LF', 'HorizontalAlignment','center','FontSize',12,'FontWeight','bold','Color',[0 0.25 0.8]); end
if any(idxH), text(mean(HF_band), ypos, 'HF', 'HorizontalAlignment','center','FontSize',12,'FontWeight','bold','Color',[0.85 0.45 0]); end

xlim([0 0.5]);
hold off;
else
    error('RR_resampled not available. Cannot perform spectral analysis.');
end