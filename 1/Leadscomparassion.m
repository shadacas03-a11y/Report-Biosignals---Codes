%% === ECG SIGNAL ANALYSIS (3 LEADS) ===
clear; close all; clc;

% === Data file sources ===
file1 = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';
%% === INTERVAL ANALYSIS: QQ, RR, SS, TT ===

% Ensure all arrays are column vectors and sorted
Q_locs = Q_locs(~isnan(Q_locs)); 
R_locs = R_locs(~isnan(R_locs)); 
S_locs = S_locs(~isnan(S_locs)); 
T_locs = T_locs(~isnan(T_locs));

% Compute time vectors
Q_times = t(Q_locs);
R_times = t(R_locs);
S_times = t(S_locs);
T_times = t(T_locs);

% Compute intervals in seconds
QQ_intervals = diff(Q_times);
RR_intervals = diff(R_times);
SS_intervals = diff(S_times);
TT_intervals = diff(T_times);

% Combine for boxplot
intervals = [QQ_intervals(:), RR_intervals(:), SS_intervals(:), TT_intervals(:)];
labels = {'QQ', 'RR', 'SS', 'TT'};

% === BOX PLOT OF INTERVALS ===
figure('Position',[100 100 800 400]);
boxplot(intervals, 'Labels', labels);
xlabel('Interval Type');
ylabel('Duration (s)');
title('Distribution of ECG Intervals (QQ, RR, SS, TT)');
grid on;

% === STATISTICS SUMMARY ===
mean_QQ = mean(QQ_intervals);
mean_RR = mean(RR_intervals);
mean_SS = mean(SS_intervals);
mean_TT = mean(TT_intervals);

fprintf('\n--- Mean Intervals ---\n');
fprintf('QQ: %.3f s\n', mean_QQ);
fprintf('RR: %.3f s\n', mean_RR);
fprintf('SS: %.3f s\n', mean_SS);
fprintf('TT: %.3f s\n', mean_TT);
fprintf('Estimated Heart Rate: %.1f bpm\n', 60/mean_RR);
