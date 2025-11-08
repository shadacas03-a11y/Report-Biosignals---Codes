# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 16:42:29 2025

@author: sydne
"""
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, find_peaks

# ---------- SETTINGS ----------
DATA_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Data"
OUTPUT_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Results/Time Window Comparisons"
FS = 100  # Sampling frequency in Hz
time_windows = [(0,30), (30,60)]  # List of time windows in seconds [(start,end), ...]

OUTPUT_DIR = os.path.join(OUTPUT_DIR, "output_windows")
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLORS = {"resting": "tab:blue", "rest": "tab:blue", "exercise": "tab:red"}

# ---------- FILTER HELPERS ----------
def butter_bandpass(lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low, high = lowcut / nyq, highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def bandpass_filter(data, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order)
    return filtfilt(b, a, data)

# ---------- OUTLIER REMOVAL ----------
def remove_outliers(data):
    if len(data) == 0:
        return data
    q1, q3 = np.percentile(data, [25,75])
    iqr = q3 - q1
    lower = q1 - 1.5*iqr
    upper = q3 + 1.5*iqr
    return data[(data >= lower) & (data <= upper)]

# ---------- GROUP FILES BY SUBJECT ----------
files = [f for f in os.listdir(DATA_DIR) if f.endswith(".txt")]
subjects = {}
for file in files:
    match = re.match(r"([AB])_", file, re.IGNORECASE)
    subj = match.group(1).upper() if match else "UNKNOWN"
    subjects.setdefault(subj, []).append(file)

# ---------- MAIN ANALYSIS ----------
for subj, subj_files in subjects.items():
    for file in subj_files:
        file_path = os.path.join(DATA_DIR, file)
        condition = "exercise" if "exercise" in file.lower() else "rest"
        color = COLORS.get(condition, "gray")

        # Load data
        data = np.loadtxt(file_path, comments="#")
        t = np.arange(len(data))/FS
        ecg = data[:,7]*1000
        resp = data[:,6]*1000
        ecg_filt = bandpass_filter(ecg, 0.5, 40, FS)
        resp_filt = bandpass_filter(resp, 0.05, 0.7, FS)

        # Prepare lists to store means per window
        mean_hr_per_window = []
        mean_br_per_window = []
        window_labels = []

        # Loop over time windows
        for (t_start, t_end) in time_windows:
            idx_start = int(t_start*FS)
            idx_end = int(t_end*FS)
            t_window = t[idx_start:idx_end]

            # ---------- ECG ----------
            ecg_win = ecg[idx_start:idx_end]
            ecg_filt_win = ecg_filt[idx_start:idx_end]
            peaks_ecg = find_peaks(ecg_filt_win, distance=FS*0.5)[0]
            rr_intervals = np.diff(t_window[peaks_ecg]) if len(peaks_ecg)>1 else np.array([np.nan])
            rr_clean = remove_outliers(rr_intervals)
            mean_hr = 60/np.mean(rr_clean) if len(rr_clean)>0 else np.nan
            mean_hr_per_window.append(mean_hr)

            # ---------- RESPIRATION ----------
            resp_win = resp[idx_start:idx_end]
            resp_filt_win = resp_filt[idx_start:idx_end]
            peaks_resp = find_peaks(resp_filt_win, distance=FS*1.0)[0]
            breath_intervals = np.diff(t_window[peaks_resp]) if len(peaks_resp)>1 else np.array([np.nan])
            breath_clean = remove_outliers(breath_intervals)
            mean_br = 60/np.mean(breath_clean) if len(breath_clean)>0 else np.nan
            mean_br_per_window.append(mean_br)

            window_labels.append(f"{t_start}-{t_end}s")

            # ---------- FIGURES ----------
            # ECG
            fig_ecg, axes_ecg = plt.subplots(1,2, figsize=(12,4))
            fig_ecg.suptitle(f"{subj} — ECG ({condition.capitalize()}, {t_start}-{t_end}s)", fontsize=14, fontweight='bold')
            axes_ecg[0].plot(t_window, ecg_win, color=color)
            axes_ecg[0].set_title("Raw ECG")
            axes_ecg[0].set_xlabel("Time [s]"); axes_ecg[0].set_ylabel("Amplitude [mV]"); axes_ecg[0].grid(True, linestyle=":", alpha=0.5)
            axes_ecg[1].plot(t_window, ecg_filt_win, color=color)
            axes_ecg[1].plot(t_window[peaks_ecg], ecg_filt_win[peaks_ecg], 'ro', label='R peaks')
            axes_ecg[1].set_title(f"Filtered ECG — Mean HR: {mean_hr:.1f} BPM")
            axes_ecg[1].set_xlabel("Time [s]"); axes_ecg[1].set_ylabel("Amplitude [mV]"); axes_ecg[1].grid(True, linestyle=":", alpha=0.5); axes_ecg[1].legend()
            fig_ecg.tight_layout(rect=[0,0,1,0.95])
            fig_ecg.savefig(os.path.join(OUTPUT_DIR, f"{subj}_ECG_{condition}_{t_start}-{t_end}s.png"), dpi=300)
            plt.close(fig_ecg)

            # Respiration
            fig_resp, axes_resp = plt.subplots(1,2, figsize=(12,4))
            fig_resp.suptitle(f"{subj} — Respiration ({condition.capitalize()}, {t_start}-{t_end}s)", fontsize=14, fontweight='bold')
            axes_resp[0].plot(t_window, resp_win, color=color)
            axes_resp[0].set_title("Raw Respiration"); axes_resp[0].set_xlabel("Time [s]"); axes_resp[0].set_ylabel("Amplitude [mV]"); axes_resp[0].grid(True, linestyle=":", alpha=0.5)
            axes_resp[1].plot(t_window, resp_filt_win, color=color)
            axes_resp[1].plot(t_window[peaks_resp], resp_filt_win[peaks_resp], 'ro', label='Peaks')
            axes_resp[1].set_title(f"Filtered Resp — Mean BR: {mean_br:.1f} BPM")
            axes_resp[1].set_xlabel("Time [s]"); axes_resp[1].set_ylabel("Amplitude [mV]"); axes_resp[1].grid(True, linestyle=":", alpha=0.5); axes_resp[1].legend()
            fig_resp.tight_layout(rect=[0,0,1,0.95])
            fig_resp.savefig(os.path.join(OUTPUT_DIR, f"{subj}_Resp_{condition}_{t_start}-{t_end}s.png"), dpi=300)
            plt.close(fig_resp)

            # Boxplots
            plt.figure(figsize=(6,4))
            plt.boxplot(rr_clean, labels=[condition])
            plt.ylabel("RR Interval (s)")
            plt.title(f"{subj} — ECG RR Interval ({t_start}-{t_end}s, outliers removed)")
            plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, f"{subj}_RR_Box_{condition}_{t_start}-{t_end}s.png"), dpi=300)
            plt.close()

            plt.figure(figsize=(6,4))
            plt.boxplot(breath_clean, labels=[condition])
            plt.ylabel("Breath Interval (s)")
            plt.title(f"{subj} — Breath Interval ({t_start}-{t_end}s, outliers removed)")
            plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, f"{subj}_Breath_Box_{condition}_{t_start}-{t_end}s.png"), dpi=300)
            plt.close()

        # ---------- SUMMARY BAR PLOTS ACROSS WINDOWS ----------
        # Heart Rate
        plt.figure(figsize=(6,4))
        plt.bar(window_labels, mean_hr_per_window, color='tab:red', alpha=0.7)
        plt.ylabel("Mean HR (BPM)")
        plt.title(f"{subj} — Mean Heart Rate Across Windows ({condition})")
        plt.grid(axis='y', linestyle=":", alpha=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, f"{subj}_MeanHR_Comparison_{condition}.png"), dpi=300)
        plt.close()

        # Breathing Rate
        plt.figure(figsize=(6,4))
        plt.bar(window_labels, mean_br_per_window, color='tab:blue', alpha=0.7)
        plt.ylabel("Mean BR (BPM)")
        plt.title(f"{subj} — Mean Breathing Rate Across Windows ({condition})")
        plt.grid(axis='y', linestyle=":", alpha=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, f"{subj}_MeanBR_Comparison_{condition}.png"), dpi=300)
        plt.close()

