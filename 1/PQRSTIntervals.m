%% === ECG ZOOM P-Q-R-S-T Detection and Interval Analysis ===
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

%% === PLOT ZOOMED SIGNAL WITH DETECTION ===
figure('Position',[100 100 1400 400]);
plot(t_zoom, ECG_zoom, 'b', 'LineWidth', 1.3); hold on;

plot(t_zoom(R_locs), ECG_zoom(R_locs), 'ro', 'MarkerFaceColor','r','MarkerSize',8);
plot(t_zoom(Q_locs), ECG_zoom(Q_locs), 'co', 'MarkerFaceColor','c','MarkerSize',8);
plot(t_zoom(S_locs), ECG_zoom(S_locs), 'mo', 'MarkerFaceColor','m','MarkerSize',8);
plot(t_zoom(P_locs), ECG_zoom(P_locs), 'go', 'MarkerFaceColor','g','MarkerSize',8);
plot(t_zoom(T_locs), ECG_zoom(T_locs), 'yo', 'MarkerFaceColor','y','MarkerSize',8);

offset = 0.3;
for i=1:length(P_locs), text(t_zoom(P_locs(i)), ECG_zoom(P_locs(i))+offset,'P','Color','g','FontWeight','bold','HorizontalAlignment','center'); end
for i=1:length(Q_locs), text(t_zoom(Q_locs(i)), ECG_zoom(Q_locs(i))-offset,'Q','Color','c','FontWeight','bold','HorizontalAlignment','center'); end
for i=1:length(R_locs), text(t_zoom(R_locs(i)), ECG_zoom(R_locs(i))+offset,'R','Color','r','FontWeight','bold','HorizontalAlignment','center'); end
for i=1:length(S_locs), text(t_zoom(S_locs(i)), ECG_zoom(S_locs(i))-offset,'S','Color','m','FontWeight','bold','HorizontalAlignment','center'); end
for i=1:length(T_locs), text(t_zoom(T_locs(i)), ECG_zoom(T_locs(i))+offset,'T','Color','y','FontWeight','bold','HorizontalAlignment','center'); end

xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('ECG Zoom: P-Q-R-S-T Detection');
legend('ECG','R Wave','Q Wave','S Wave','P Wave','T Wave','Location','best');
grid on;
xlim([zoom_start zoom_end]);
ylim([min(ECG_zoom)-0.5 max(ECG_zoom)+0.5]);

disp('✅ ECG zoom detection (P-Q-R-S-T) completed successfully.');

%% === CALCULATE INTERVALS AND ENHANCED BOXPLOTS (SAFE VERSION) ===
% filtered intervals calculation of NaN
RR_intervals = diff(t_zoom(R_locs(~isnan(R_locs))));
QQ_intervals = diff(t_zoom(Q_locs(~isnan(Q_locs))));
SS_intervals = diff(t_zoom(S_locs(~isnan(S_locs))));
TT_intervals = diff(t_zoom(T_locs(~isnan(T_locs))));

% --- Filled with NaN---
max_len = max([length(RR_intervals), length(QQ_intervals), length(SS_intervals), length(TT_intervals)]);
RR_pad = [RR_intervals, NaN(1,max_len-length(RR_intervals))];
QQ_pad = [QQ_intervals, NaN(1,max_len-length(QQ_intervals))];
SS_pad = [SS_intervals, NaN(1,max_len-length(SS_intervals))];
TT_pad = [TT_intervals, NaN(1,max_len-length(TT_intervals))];

% Concatenaded matrix
intervals_matrix = [RR_pad', QQ_pad', SS_pad', TT_pad'];

%% === ENHANCED BOXPLOT ===
figure('Position',[200 200 1000 500]);
boxplot(intervals_matrix, {'RR','QQ','SS','TT'}, 'Whisker',1.5, 'Colors','k', 'Symbol','r+');
hold on;

% Means and median calculation
means = mean(intervals_matrix,1,'omitnan');
medians = median(intervals_matrix,1,'omitnan');

% Plot median and mean
for i = 1:size(intervals_matrix,2)
    plot(i, means(i), 'b*', 'MarkerSize',10, 'LineWidth',1.5);  % media
    plot(i, medians(i), 'gs', 'MarkerSize',8, 'LineWidth',1.5); % mediana
end

ylabel('Interval Duration (s)');
title('Enhanced Boxplots of ECG Intervals (RR, QQ, SS, TT)');
legend({'Outliers','Mean','Median'}, 'Location','best');
grid on;

disp('✅ ECG intervals calculated and enhanced boxplots generated successfully.');
