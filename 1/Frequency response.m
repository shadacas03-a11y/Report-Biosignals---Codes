% === Data file source ===
filename = 'C:\Users\CAMILA CASTRO\Downloads\Application to cardiovascular fitness\1. Sharon_LEAD1_resting_converted.txt';

% Open the file
fid = fopen(filename, 'r');
if fid == -1
    error('No se pudo abrir el archivo. Verifica la ruta.');
end

% Skip header lines until "# EndOfHeader"
while true
    line = fgetl(fid);
    if ~ischar(line)
        error('No se encontró "# EndOfHeader" en el archivo.');
    end
    if contains(line, '# EndOfHeader')
        break;
    end
end

% Read the numerical data
data = textscan(fid, '%f%f%f%f%f%f', 'Delimiter', '\t', 'CollectOutput', true);
fclose(fid);

data = data{1};  
fs = 100;  % Sampling frequency [Hz]

% Select ECG channel A2 (Lead II)
ECG = data(:,6);
t = (0:length(ECG)-1)/fs;  % time vector



%%==== Frequency response ==== %%
% High-pass filter (<0.5 Hz)
[b_hp, a_hp] = butter(2, 0.5/(fs/2), 'high');

% Low-pass filter (>40 Hz)
[b_lp, a_lp] = butter(2, 40/(fs/2), 'low');
% High-pass frequency response
[H_hp, f_hp] = freqz(b_hp, a_hp, 1024, fs);
figure;
plot(f_hp, abs(H_hp), 'r', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('High-pass filter frequency response (<0.5 Hz)');
grid on;

% Low-pass frequency response
[H_lp, f_lp] = freqz(b_lp, a_lp, 1024, fs);
figure;
plot(f_lp, abs(H_lp), 'b', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Low-pass filter frequency response (<40 Hz)');
grid on;

