# -*- coding: utf-8 -*-
"""
Created on Fri Nov 7 2025

ECG and Pulse Sensor RR-interval analysis for 30-second windows.
Computes R peaks, RR intervals, HR, saves diagnostic plots with separate subplots,
and prepares data for Bland–Altman plots.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, find_peaks

# ---------- SETTINGS ----------
DATA_DIR = r"C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Data"
OUTPUT_DIR = r"C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Results/ECG_Pulse_BA"
FS = 100  # Hz
WINDOW_DURATION = 30  # seconds
window_samples = FS * WINDOW_DURATION

COLORS = {
    "rest": "tab:blue",
    "exercise": "tab:red"
}

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------- FILTER HELPERS ----------
def bandpass_filter(data, lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    b, a = butter(order, [lowcut/nyq, highcut/nyq], btype='band')
    return filtfilt(b, a, data)

# ---------- HR & RR FUNCTION ----------
def compute_hr(signal_window, fs):
    peaks, _ = find_peaks(signal_window, distance=int(0.5*fs))
    rr_intervals = np.diff(peaks)/fs if len(peaks) > 1 else np.array([np.nan])
    mean_hr = 60/np.mean(rr_intervals) if len(rr_intervals) > 0 else np.nan
    return mean_hr, rr_intervals, peaks

# ---------- MAIN LOOP ----------
files = [f for f in os.listdir(DATA_DIR) if f.endswith(".txt")]
summary_records = []

for file in files:
    print(f"Processing {file}...")
    
    # Subject and condition
    subj_match = re.search(r"([AB])_", file)
    cond_match = re.search(r"(rest|exercise)", file.lower())
    subject = subj_match.group(1) if subj_match else "Unknown"
    condition = cond_match.group(1) if cond_match else "unspecified"
    color = COLORS.get(condition, "gray")
    
    # Load raw data
    df = np.loadtxt(os.path.join(DATA_DIR, file), comments="#")
    ecg = df[:, 7] * 1000        # ECG in mV 
    pulse = (df[:, 5]  / 1023.0) * 3.3 * 1000
  
    
    # Filter signals once
    ecg_filt = bandpass_filter(ecg, 0.5, 40, FS)
    if pulse is not None:
        pulse_filt = bandpass_filter(pulse, 0.5, 10, FS)
    
    # Windowing
    total_samples = len(ecg)
    num_windows = int(np.ceil(total_samples/window_samples))
    
    for i in range(num_windows):
        start = i*window_samples
        end = min((i+1)*window_samples, total_samples)
        t = np.arange(end-start)/FS
        
        ecg_win = ecg_filt[start:end]
        if pulse is not None:
            pulse_win = pulse_filt[start:end]
        
        # Compute HR
        hr_ecg, rr_ecg, peaks_ecg = compute_hr(ecg_win, FS)
        if pulse is not None:
            hr_pulse, rr_pulse, peaks_pulse = compute_hr(pulse_win, FS)
        else:
            hr_pulse = np.nan
            peaks_pulse = np.array([])
        
        # Save summary
        summary_records.append({
            'Subject': subject,
            'Condition': condition,
            'Window': i+1,
            'Start_Time_s': start/FS,
            'End_Time_s': end/FS,
            'HR_ECG': round(hr_ecg,2),
            'HR_Pulse': round(hr_pulse,2),
            'Num_R_Peaks_ECG': len(peaks_ecg),
            'Num_R_Peaks_Pulse': len(peaks_pulse)
        })
        
        # ---------- Diagnostic plot with separate subplots ----------
        plt.figure(figsize=(12,6))
        
        # ECG subplot
        plt.subplot(2,1,1)
        plt.plot(t, ecg_win, 'k', label='ECG')
        plt.plot(t[peaks_ecg], ecg_win[peaks_ecg], 'ro', label='R peaks')
        plt.title(f"{subject} {condition} — Window {i+1} — ECG HR: {hr_ecg:.1f} BPM")
        plt.ylabel("ECG [mV]")
        plt.legend()
        plt.grid(True, linestyle=":")
        
        # Pulse subplot
        plt.subplot(2,1,2)
        if pulse is not None:
            plt.plot(t, pulse_win, 'b', label='Pulse')
            plt.plot(t[peaks_pulse], pulse_win[peaks_pulse], 'go', label='Pulse peaks')
            plt.title(f"{subject} {condition} — Pulse HR: {hr_pulse:.1f} BPM")
        else:
            plt.plot([], [], 'b', label='Pulse not available')
            plt.title(f"{subject} {condition} — Pulse data missing")
        plt.xlabel("Time [s]")
        plt.ylabel("Amplitude")
        plt.legend()
        plt.grid(True, linestyle=":")
        
        plt.tight_layout()
        plot_file = os.path.join(OUTPUT_DIR, f"{subject}_{condition}_win{i+1}.png")
        plt.savefig(plot_file, dpi=300)
        plt.close()

# Save summary CSV
summary_df = pd.DataFrame(summary_records)
csv_file = os.path.join(OUTPUT_DIR, "ECG_Pulse_HR_summary.csv")
summary_df.to_csv(csv_file, index=False)
print(f"\n✅ Summary saved to: {csv_file}")

# ---------- Bland–Altman Plot ----------
# Only include windows where pulse data exists
ba_df = summary_df.dropna(subset=['HR_Pulse']).copy()
ba_df['Mean_HR'] = (ba_df['HR_ECG'] + ba_df['HR_Pulse'])/2
ba_df['Diff_HR'] = ba_df['HR_ECG'] - ba_df['HR_Pulse']

mean_diff = np.mean(ba_df['Diff_HR'])
sd_diff = np.std(ba_df['Diff_HR'])

plt.figure(figsize=(8,6))
plt.scatter(ba_df['Mean_HR'], ba_df['Diff_HR'], color='tab:purple')
plt.axhline(mean_diff, color='red', label='Mean Diff')
plt.axhline(mean_diff + 1.96*sd_diff, color='gray', linestyle='--', label='95% Limits')
plt.axhline(mean_diff - 1.96*sd_diff, color='gray', linestyle='--')
plt.xlabel('Mean HR (BPM)')
plt.ylabel('Difference HR (ECG - Pulse)')
plt.title('Bland–Altman Plot')
plt.legend()
plt.grid(True, linestyle=':')
plt.tight_layout()
ba_plot_file = os.path.join(OUTPUT_DIR, "BlandAltman_HR.png")
plt.savefig(ba_plot_file, dpi=300)
plt.show()
