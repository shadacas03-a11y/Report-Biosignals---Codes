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

%% === ZOOM INTERVAL (Optional) ===
zoom_start = 176; 
zoom_end = 182;
idx_zoom = t >= zoom_start & t <= zoom_end;
t_zoom = t(idx_zoom);
ECG_zoom = ECG_smooth(idx_zoom);

%% === P-Q-R-S-T DETECTION ===
% Find R peaks
minPeakHeight = 0.25 * max(ECG_zoom);
minPeakDistance = 0.5 * fs;
[R_peaks, R_locs] = findpeaks(ECG_zoom,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);

% Find minima (for Q and S)
[~, min_locs] = findpeaks(-ECG_zoom);

P_locs=[]; Q_locs=[]; S_locs=[]; T_locs=[];

for i=1:length(R_locs)
    r = R_locs(i);

    % Q-wave
    q_candidates = min_locs(min_locs<r);
    if ~isempty(q_candidates)
        q = q_candidates(end);
        Q_locs = [Q_locs q];
    else
        q = NaN;
    end

    % S-wave
    s_candidates = min_locs(min_locs>r);
    if ~isempty(s_candidates)
        s = s_candidates(1);
        S_locs = [S_locs s];
    else
        s = NaN;
    end

    % P-wave
    if ~isnan(q)
        region = max(q-round(0.25*fs),1):q;
        [pksP, locsP] = findpeaks(ECG_zoom(region));
        if ~isempty(pksP)
            P_locs = [P_locs region(locsP(end))];
        end
    end

    % T-wave
    if ~isnan(s)
        region = s:min(s+round(0.5*fs), length(ECG_zoom));
        [pksT, locsT] = findpeaks(ECG_zoom(region));
        if ~isempty(pksT)
            T_locs = [T_locs region(locsT(1))];
        end
    end
end

%% === INTERVAL CALCULATION ===
RR_intervals = diff(t_zoom(R_locs(~isnan(R_locs))));
QQ_intervals = diff(t_zoom(Q_locs(~isnan(Q_locs))));
SS_intervals = diff(t_zoom(S_locs(~isnan(S_locs))));
TT_intervals = diff(t_zoom(T_locs(~isnan(T_locs))));

%% === REMOVE OUTLIERS ===
remove_outliers = @(data) data(data >= prctile(data,25)-1.5*(prctile(data,75)-prctile(data,25)) & ...
                                data <= prctile(data,75)+1.5*(prctile(data,75)-prctile(data,25)));

RR_clean = remove_outliers(RR_intervals);
QQ_clean = remove_outliers(QQ_intervals);
SS_clean = remove_outliers(SS_intervals);
TT_clean = remove_outliers(TT_intervals);

%% === STATISTICS AND NORMALITY TEST ===
intervals_data = {RR_clean, QQ_clean, SS_clean, TT_clean};
names = {'RR','QQ','SS','TT'};

fprintf('\n--- Interval Statistics ---\n');
for i = 1:length(intervals_data)
    data = intervals_data{i};
    if isempty(data)
        fprintf('%s-interval: No valid data.\n', names{i});
        continue;
    end
    fprintf('%s: mean = %.3f s, std = %.3f s, median = %.3f s\n', ...
        names{i}, mean(data), std(data), median(data));
end

fprintf('\n--- Normality Test (Lilliefors) ---\n');
for i = 1:length(intervals_data)
    data = intervals_data{i};
    if length(data) < 4
        fprintf('%s: Not enough data (N = %d)\n', names{i}, length(data));
        continue;
    end
    [h,p] = lillietest(data);
    if h==0
        fprintf('%s: p = %.3f → normal distribution\n', names{i}, p);
    else
        fprintf('%s: p = %.3f → NOT normal\n', names{i}, p);
    end
end

%% === RESAMPLE RR INTERVALS AT 4 Hz ===
if ~isempty(RR_clean)
    t_RR = cumsum(RR_clean);
    fs_resample = 4;
    t_uniform = t_RR(1):1/fs_resample:t_RR(end);
    RR_resampled = interp1(t_RR, RR_clean, t_uniform, 'pchip');

    figure; 
    plot(t_uniform, RR_resampled, '-o','LineWidth',1.2,'MarkerSize',4);
    xlabel('Time (s)'); ylabel('RR Interval (s)');
    title('RR Intervals Resampled at 4 Hz'); grid on;
else
    warning('RR_clean unavailable');
end

%% === SPECTRAL ANALYSIS ===
if exist('RR_resampled','var') && ~isempty(RR_resampled)
    RR_detrended = detrend(RR_resampled);
    [Pxx,f] = pwelch(RR_detrended,[],[],[],fs_resample);

    LF_band = [0.04 0.15];
    HF_band = [0.15 0.4];

    LF_power = bandpower(Pxx,f,LF_band,'psd');
    HF_power = bandpower(Pxx,f,HF_band,'psd');
    LF_HF_ratio = LF_power / HF_power;

    fprintf('\nSpectral Analysis:\n');
    fprintf('LF Power = %.4f s^2\n', LF_power);
    fprintf('HF Power = %.4f s^2\n', HF_power);
    fprintf('LF/HF Ratio = %.2f\n', LF_HF_ratio);

    figure;
    plot(f,Pxx,'LineWidth',1.5); hold on;
    xline(LF_band,'r--','LF band'); xline(HF_band,'g--','HF band');
    xlabel('Frequency (Hz)'); ylabel('Power Spectral Density');
    title('PSD of RR Intervals (4 Hz Resampled)'); grid on;
end
