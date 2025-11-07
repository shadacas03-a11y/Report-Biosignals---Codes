%% === ECG ZOOM P-Q-R-S-T Detection and Interval Analysis (Complete Script) ===
clc; clear; close all;

%% === LOAD DATA ===
filename = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file. Check the path.');
end

% Skip header lines until '# EndOfHeader'
while true
    line = fgetl(fid);
    if ~ischar(line)
        error('# EndOfHeader not found.');
    end
    if contains(line, '# EndOfHeader')
        break;
    end
end

data = textscan(fid, '%f%f%f%f%f%f', 'Delimiter', '\t', 'CollectOutput', true);
fclose(fid);
data = data{1};

fs = 100;                     % Sampling frequency [Hz]
ECG = data(:,6);              % Lead II
t = (0:length(ECG)-1)/fs;

%% === FILTERING ===
[b_hp,a_hp] = butter(2,0.5/(fs/2),'high');
[b_lp,a_lp] = butter(4,35/(fs/2),'low');
ECG_filt = filtfilt(b_lp,a_lp,filtfilt(b_hp,a_hp,ECG))*1000;  % mV
ECG_smooth = smoothdata(ECG_filt,'gaussian',round(0.015*fs));

%% === ZOOM INTERVAL ===
zoom_start = 176;  % seconds
zoom_end   = 182;  % seconds
idx_zoom = t >= zoom_start & t <= zoom_end;
t_zoom = t(idx_zoom);
ECG_zoom = ECG_smooth(idx_zoom);

%% === DERIVATIVE-BASED P-Q-R-S-T DETECTION ===
% Find R peaks
minPeakHeight = 0.25 * max(ECG_zoom);
minPeakDistance = 0.5 * fs;
[R_peaks, R_locs] = findpeaks(ECG_zoom, 'MinPeakHeight',minPeakHeight, 'MinPeakDistance',minPeakDistance);

% Find all minima for Q and S
[~, min_locs] = findpeaks(-ECG_zoom);

% Initialize arrays
P_locs = []; Q_locs = []; S_locs = []; T_locs = [];

for i = 1:length(R_locs)
    r = R_locs(i);

    % --- Q: nearest minimum before R ---
    q_candidates = min_locs(min_locs < r);
    if ~isempty(q_candidates)
        q = q_candidates(end);
        Q_locs = [Q_locs q];
    else
        q = NaN;
    end

    % --- S: nearest minimum after R ---
    s_candidates = min_locs(min_locs > r);
    if ~isempty(s_candidates)
        s = s_candidates(1);
        S_locs = [S_locs s];
    else
        s = NaN;
    end

    % --- P: local maximum before Q ---
    if ~isnan(q)
        search_region = max(q - round(0.25*fs), 1):q;
        [pksP, locsP] = findpeaks(ECG_zoom(search_region));
        if ~isempty(pksP)
            P_locs = [P_locs search_region(locsP(end))];
        end
    end

    % --- T: local maximum after S ---
    if ~isnan(s)
        search_region = s:min(s + round(0.5*fs), length(ECG_zoom));
        [pksT, locsT] = findpeaks(ECG_zoom(search_region));
        if ~isempty(pksT)
            T_locs = [T_locs search_region(locsT(1))];
        end
    end
end

%% === CALCULATE INTERVALS ===
RR_intervals = diff(t_zoom(R_locs(~isnan(R_locs))));
QQ_intervals = diff(t_zoom(Q_locs(~isnan(Q_locs))));
SS_intervals = diff(t_zoom(S_locs(~isnan(S_locs))));
TT_intervals = diff(t_zoom(T_locs(~isnan(T_locs))));

%% === ECG INTERVALS ANALYSIS WITH BOXPLOTS AND OUTLIER REMOVAL ===
% This section detects P-Q-R-S-T waves in the ECG, calculates RR, QQ, SS, TT intervals,
% and visualizes their distributions using boxplots before and after removing outliers.

% --- Original intervals ---
all_intervals = [RR_intervals, QQ_intervals, SS_intervals, TT_intervals];
group = [repmat({'RR'},1,length(RR_intervals)), ...
         repmat({'QQ'},1,length(QQ_intervals)), ...
         repmat({'SS'},1,length(SS_intervals)), ...
         repmat({'TT'},1,length(TT_intervals))];

figure('Position',[100 100 1400 500]);

% --- Boxplot BEFORE outlier removal ---
subplot(1,2,1);
boxplot(all_intervals, group, 'Whisker',1.5, 'Colors','k');
hold on;

% Scatter all points for clarity
scatter(ones(size(RR_intervals)), RR_intervals, 50, 'r', 'filled');
scatter(2*ones(size(QQ_intervals)), QQ_intervals, 50, 'r', 'filled');
scatter(3*ones(size(SS_intervals)), SS_intervals, 50, 'r', 'filled');
scatter(4*ones(size(TT_intervals)), TT_intervals, 50, 'r', 'filled');

% Mark mean (blue *) and median (green square)
plot(1, mean(RR_intervals), 'b*', 'MarkerSize',10);
plot(2, mean(QQ_intervals), 'b*', 'MarkerSize',10);
plot(3, mean(SS_intervals), 'b*', 'MarkerSize',10);
plot(4, mean(TT_intervals), 'b*', 'MarkerSize',10);

plot(1, median(RR_intervals), 'gs', 'MarkerSize',8);
plot(2, median(QQ_intervals), 'gs', 'MarkerSize',8);
plot(3, median(SS_intervals), 'gs', 'MarkerSize',8);
plot(4, median(TT_intervals), 'gs', 'MarkerSize',8);

title('Boxplots Before Outlier Removal');
ylabel('Interval Duration (s)');
grid on;

% --- Remove outliers using 1.5*IQR rule ---
remove_outliers = @(data) data(data >= prctile(data,25)-1.5*(prctile(data,75)-prctile(data,25)) & ...
                                data <= prctile(data,75)+1.5*(prctile(data,75)-prctile(data,25)));

RR_clean = remove_outliers(RR_intervals);
QQ_clean = remove_outliers(QQ_intervals);
SS_clean = remove_outliers(SS_intervals);
TT_clean = remove_outliers(TT_intervals);

disp(['RR intervals: ', num2str(length(RR_intervals)), ' -> ', num2str(length(RR_clean))]);
disp(['QQ intervals: ', num2str(length(QQ_intervals)), ' -> ', num2str(length(QQ_clean))]);
disp(['SS intervals: ', num2str(length(SS_intervals)), ' -> ', num2str(length(SS_clean))]);
disp(['TT intervals: ', num2str(length(TT_intervals)), ' -> ', num2str(length(TT_clean))]);

all_intervals_clean = [RR_clean, QQ_clean, SS_clean, TT_clean];
group_clean = [repmat({'RR'},1,length(RR_clean)), ...
               repmat({'QQ'},1,length(QQ_clean)), ...
               repmat({'SS'},1,length(SS_clean)), ...
               repmat({'TT'},1,length(TT_clean))];

% --- Boxplot AFTER outlier removal ---
subplot(1,2,2);
boxplot(all_intervals_clean, group_clean, 'Whisker',1.5, 'Colors','k');
hold on;

% Scatter cleaned points
scatter(ones(size(RR_clean)), RR_clean, 50, 'r', 'filled');
scatter(2*ones(size(QQ_clean)), QQ_clean, 50, 'r', 'filled');
scatter(3*ones(size(SS_clean)), SS_clean, 50, 'r', 'filled');
scatter(4*ones(size(TT_clean)), TT_clean, 50, 'r', 'filled');

% Mean & median markers
plot(1, mean(RR_clean), 'b*', 'MarkerSize',10);
plot(2, mean(QQ_clean), 'b*', 'MarkerSize',10);
plot(3, mean(SS_clean), 'b*', 'MarkerSize',10);
plot(4, mean(TT_clean), 'b*', 'MarkerSize',10);

plot(1, median(RR_clean), 'gs', 'MarkerSize',8);
plot(2, median(QQ_clean), 'gs', 'MarkerSize',8);
plot(3, median(SS_clean), 'gs', 'MarkerSize',8);
plot(4, median(TT_clean), 'gs', 'MarkerSize',8);

title('Boxplots After Outlier Removal');
ylabel('Interval Duration (s)');
grid on;

sgtitle('ECG Intervals: Comparison Before and After Outlier Removal');

disp('✅ Boxplots with scatter and outlier removal completed successfully.');


