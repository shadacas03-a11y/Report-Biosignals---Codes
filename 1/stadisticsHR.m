%% === ECG P-Q-R-S-T Detection, Interval Analysis, and Heart Rate Evolution ===
clc; clear; close all;

%% === LOAD DATA ===
filename = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';
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

fs = 100;                   % Sampling frequency [Hz]
ECG = data(:,6);            % Lead II
t = (0:length(ECG)-1)/fs;

%% === FILTERING ===
[b_hp,a_hp] = butter(2,0.5/(fs/2),'high');
[b_lp,a_lp] = butter(4,35/(fs/2),'low');
ECG_filt = filtfilt(b_lp,a_lp,filtfilt(b_hp,a_hp,ECG))*1000;  % mV
ECG_smooth = smoothdata(ECG_filt,'gaussian',round(0.015*fs));

%% === ZOOM INTERVAL (Optional, small segment) ===
zoom_start = 176; zoom_end = 182;
idx_zoom = t >= zoom_start & t <= zoom_end;
t_zoom = t(idx_zoom);
ECG_zoom = ECG_smooth(idx_zoom);

%% === P-Q-R-S-T DETECTION (Derivative-Based) ===
% Find R peaks
minPeakHeight = 0.25 * max(ECG_zoom);
minPeakDistance = 0.5 * fs;
[R_peaks, R_locs] = findpeaks(ECG_zoom,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);

% Find all minima for Q and S
[~, min_locs] = findpeaks(-ECG_zoom);

% Initialize arrays
P_locs=[]; Q_locs=[]; S_locs=[]; T_locs=[];

for i=1:length(R_locs)
    r = R_locs(i);

    % Q: nearest minimum before R
    q_candidates = min_locs(min_locs<r);
    if ~isempty(q_candidates)
        q = q_candidates(end);
        Q_locs = [Q_locs q];
    else
        q = NaN;
    end

    % S: nearest minimum after R
    s_candidates = min_locs(min_locs>r);
    if ~isempty(s_candidates)
        s = s_candidates(1);
        S_locs = [S_locs s];
    else
        s = NaN;
    end

    % P: local max before Q
    if ~isnan(q)
        search_region = max(q-round(0.25*fs),1):q;
        [pksP, locsP] = findpeaks(ECG_zoom(search_region));
        if ~isempty(pksP)
            P_locs = [P_locs search_region(locsP(end))];
        end
    end

    % T: local max after S
    if ~isnan(s)
        search_region = s:min(s+round(0.5*fs), length(ECG_zoom));
        [pksT, locsT] = findpeaks(ECG_zoom(search_region));
        if ~isempty(pksT)
            T_locs = [T_locs search_region(locsT(1))];
        end
    end
end

%% === INTERVAL CALCULATION ===
RR_intervals = diff(t_zoom(R_locs(~isnan(R_locs))));
QQ_intervals = diff(t_zoom(Q_locs(~isnan(Q_locs))));
SS_intervals = diff(t_zoom(S_locs(~isnan(S_locs))));
TT_intervals = diff(t_zoom(T_locs(~isnan(T_locs))));

%% === REMOVE OUTLIERS (Robust) ===
remove_outliers = @(data) data(data >= prctile(data,25)-1.5*(prctile(data,75)-prctile(data,25)) & ...
                                data <= prctile(data,75)+1.5*(prctile(data,75)-prctile(data,25)));

RR_clean = remove_outliers(RR_intervals);
QQ_clean = remove_outliers(QQ_intervals);
SS_clean = remove_outliers(SS_intervals);
TT_clean = remove_outliers(TT_intervals);

%% === STATISTICS AND NORMALITY TEST (Robust) ===
intervals_data = {RR_clean, QQ_clean, SS_clean, TT_clean};
names = {'RR','QQ','SS','TT'};

fprintf('\n--- Interval Statistics ---\n');
for i = 1:length(intervals_data)
    data = intervals_data{i};
    if isempty(data)
        fprintf('%s-interval: No valid data available.\n', names{i});
        continue;
    end
    fprintf('%s-interval: mean = %.3f s, std = %.3f s, median = %.3f s\n', ...
        names{i}, mean(data), std(data), median(data));
end

fprintf('\n--- Normality Test (Lilliefors) ---\n');
for i = 1:length(intervals_data)
    data = intervals_data{i};
    if length(data) < 4
        fprintf('%s-interval: Not enough data for normality test (N = %d)\n', names{i}, length(data));
        continue
    end
    [h,p] = lillietest(data);
    if h==0
        fprintf('%s-interval: p = %.3f -> data seems normal\n', names{i}, p);
    else
        fprintf('%s-interval: p = %.3f -> data does NOT seem normal\n', names{i}, p);
    end
end

%% === HEART RATE AND HRV ANALYSIS ===
if ~isempty(RR_clean)
    % --- Heart Rate series ---
    HR_series = 60 ./ RR_clean;       % bpm
    t_RR = cumsum(RR_clean);          % cumulative time

    % --- Plot HR evolution ---
    figure('Position',[100 100 800 400]);
    plot(t_RR, HR_series,'-o','LineWidth',1.2,'MarkerSize',6);
    xlabel('Time (s)');
    ylabel('Heart Rate (bpm)');
    title('Heart Rate Evolution Over Time');
    grid on;

    % --- Average HR ---
    avg_HR = mean(HR_series);
    fprintf('\nAverage Heart Rate = %.1f bpm\n', avg_HR);

    % --- SDNN (overall HR variability) ---
    SDNN = std(RR_clean);  % in seconds
    fprintf('SDNN (Standard Deviation of RR intervals) = %.4f s\n', SDNN);

    % --- HRV (RMSSD: root mean square of successive differences) ---
    if length(RR_clean) > 1
        diff_RR = diff(RR_clean);
        RMSSD = sqrt(mean(diff_RR.^2)) * 1000;  % convert to ms
        fprintf('HRV (RMSSD) = %.2f ms\n', RMSSD);
    else
        warning('Not enough RR intervals to compute HRV.');
    end

    % --- Optional: histogram of RR intervals ---
    figure('Position',[100 100 600 400]);
    histogram(RR_clean,'BinWidth',0.01,'FaceColor','b');
    xlabel('RR Interval (s)');
    ylabel('Count');
    title('Distribution of RR Intervals (Cleaned)');
    grid on;

else
    warning('RR_clean data not available. Cannot compute HR or HRV.');
end
