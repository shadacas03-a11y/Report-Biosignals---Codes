%% === LOAD DATA ===
filename = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';

fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file. Check the path.');
end

% Skip header
while true
    line = fgetl(fid);
    if ~ischar(line)
        error('# EndOfHeader not found.');
    end
    if contains(line, '# EndOfHeader')
        break;
    end
end

% Read numerical data
data = textscan(fid, '%f%f%f%f%f%f', 'Delimiter', '\t', 'CollectOutput', true);
fclose(fid);

data = data{1};
fs = 100; % sampling frequency [Hz]

% Lead II
ECG = data(:,6);
t = (0:length(ECG)-1)/fs;

%% === FILTERING ===
[b_hp, a_hp] = butter(2, 0.5/(fs/2), 'high');
ECG_hp = filtfilt(b_hp, a_hp, ECG);

[b_lp, a_lp] = butter(2, 40/(fs/2), 'low');
ECG_filtered = filtfilt(b_lp, a_lp, ECG_hp);

ECG_filtered_mV = ECG_filtered*1000;

%% === QRS DETECTION ===
minPeakHeight = 0.5 * max(ECG_filtered_mV);
minPeakDistance = 0.6 * fs;
[~, R_locs] = findpeaks(ECG_filtered_mV, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDistance);

%% === Q AND S WAVES DETECTION ===
Q_locs = [];
S_locs = [];
window_samples = round(0.05*fs); % 50 ms before and after R

for i = 1:length(R_locs)
    % Q wave: local minimum before R
    win_start = max(R_locs(i)-window_samples, 1);
    win_end = R_locs(i)-1;
    [~, q_idx] = min(ECG_filtered_mV(win_start:win_end));
    Q_locs = [Q_locs; win_start-1 + q_idx];
    
    % S wave: local minimum after R
    win_start = R_locs(i)+1;
    win_end = min(R_locs(i)+window_samples, length(ECG_filtered_mV));
    [~, s_idx] = min(ECG_filtered_mV(win_start:win_end));
    S_locs = [S_locs; win_start-1 + s_idx];
end

%% === P AND T WAVES DETECTION ===
P_locs = [];
T_locs = [];
window_samples_PT = round(0.2 * fs); % 200 ms

for i = 1:length(R_locs)
    % P wave before R
    if R_locs(i)-window_samples_PT > 0
        [~, p_idx] = max(ECG_filtered_mV(R_locs(i)-window_samples_PT:R_locs(i)-10));
        P_locs = [P_locs; R_locs(i)-window_samples_PT-1 + p_idx];
    end
    % T wave after S
    if S_locs(i)+10+window_samples_PT <= length(ECG_filtered_mV)
        [~, t_idx] = max(ECG_filtered_mV(S_locs(i)+10:S_locs(i)+window_samples_PT));
        T_locs = [T_locs; S_locs(i)+10-1 + t_idx];
    end
end

%% === ZOOM WINDOW (2.5 to 3 min) ===
zoom_start = 2.5*60;
zoom_end   = 3*60;

idx_zoom = t >= zoom_start & t <= zoom_end;

% Keep only peaks inside zoom window
R_zoom = R_locs(R_locs/fs >= 2.5 & R_locs/fs <= 3);
Q_zoom = Q_locs(Q_locs/fs >= 2.5 & Q_locs/fs <= 3);
S_zoom = S_locs(S_locs/fs >= 2.5 & S_locs/fs <= 3);
P_zoom = P_locs(P_locs/fs >= 2.5 & P_locs/fs <= 3);
T_zoom = T_locs(T_locs/fs >= 2.5 & T_locs/fs <= 3);

%% === PLOT ZOOMED ECG WITH ALL WAVES ===
figure('Position',[100 100 1400 400]);
plot(t(idx_zoom), ECG_filtered_mV(idx_zoom), 'b','LineWidth',1.2); hold on;

% Plot each wave
plot(t(P_zoom), ECG_filtered_mV(P_zoom), 'go','MarkerFaceColor','g','MarkerSize',8);
plot(t(Q_zoom), ECG_filtered_mV(Q_zoom), 'co','MarkerFaceColor','c','MarkerSize',8);
plot(t(R_zoom), ECG_filtered_mV(R_zoom), 'ro','MarkerFaceColor','r','MarkerSize',8);
plot(t(S_zoom), ECG_filtered_mV(S_zoom), 'mo','MarkerFaceColor','m','MarkerSize',8);
plot(t(T_zoom), ECG_filtered_mV(T_zoom), 'yo','MarkerFaceColor','y','MarkerSize',8);

% Add labels
for i=1:length(P_zoom)
    text(t(P_zoom(i)), ECG_filtered_mV(P_zoom(i))+0.3,'P','Color','g','FontWeight','bold','HorizontalAlignment','center');
end
for i=1:length(Q_zoom)
    text(t(Q_zoom(i)), ECG_filtered_mV(Q_zoom(i))-0.3,'Q','Color','c','FontWeight','bold','HorizontalAlignment','center');
end
for i=1:length(R_zoom)
    text(t(R_zoom(i)), ECG_filtered_mV(R_zoom(i))+0.3,'R','Color','r','FontWeight','bold','HorizontalAlignment','center');
end
for i=1:length(S_zoom)
    text(t(S_zoom(i)), ECG_filtered_mV(S_zoom(i))-0.3,'S','Color','m','FontWeight','bold','HorizontalAlignment','center');
end
for i=1:length(T_zoom)
    text(t(T_zoom(i)), ECG_filtered_mV(T_zoom(i))+0.3,'T','Color','y','FontWeight','bold','HorizontalAlignment','center');
end

xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('Zoomed ECG (Lead II) 2.5-3 min with P, Q, R, S, T waves');
grid on;
xlim([zoom_start zoom_end]);
ylim([min(ECG_filtered_mV(idx_zoom))-0.5, max(ECG_filtered_mV(idx_zoom))+0.5]);
legend('ECG','P wave','Q wave','R peak','S wave','T wave','Location','best');
