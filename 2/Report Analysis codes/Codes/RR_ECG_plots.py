# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 15:50:26 2025

@author: sydne
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, find_peaks

# ---------- SETTINGS ----------
DATA_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Data"
OUTPUT_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Results/RR Analysis"
FS = 100  # Hz
PLOT_DURATION = 5  # seconds
OUTPUT_DIR = os.path.join(OUTPUT_DIR, f"output_{PLOT_DURATION}s")
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLORS = {
    "resting": "tab:blue",
    "rest": "tab:blue",
    "exercise": "tab:red",
}

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------- FILTER HELPERS ----------
def butter_bandpass(lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low, high = lowcut / nyq, highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def bandpass_filter(data, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order)
    return filtfilt(b, a, data)

# ---------- GROUP FILES BY SUBJECT ----------
files = [f for f in os.listdir(DATA_DIR) if f.endswith(".txt")]
subjects = {}
for file in files:
    match = re.match(r"([AB])_", file, re.IGNORECASE)
    subj = match.group(1).upper() if match else "UNKNOWN"
    subjects.setdefault(subj, []).append(file)

# ---------- PROCESS EACH SUBJECT ----------
for subj, subj_files in subjects.items():

    # ----- ECG FIGURE -----
    plt.figure(figsize=(12, 8))
    plt.suptitle(f"Subject {subj} — ECG Signals with R Peaks", fontsize=14, fontweight='bold')

    for file in subj_files:
        file_path = os.path.join(DATA_DIR, file)
        data = np.loadtxt(file_path, comments="#")
        ecg = data[:, 7] * 1000  # convert to mV
        t = np.arange(len(ecg)) / FS

        condition = "exercise" if "exercise" in file.lower() else "rest"
        color = COLORS.get(condition, "gray")
        ecg_filt = bandpass_filter(ecg, 0.5, 40, FS)

        # Limit to first 30 seconds
        num_samples = min(len(ecg), FS * PLOT_DURATION)
        t_plot = t[:num_samples]
        ecg_plot = ecg[:num_samples]
        ecg_filt_plot = ecg_filt[:num_samples]

        # Detect R peaks in filtered signal
        peaks, _ = find_peaks(ecg_filt_plot, distance=FS*0.5)
        rr_intervals = np.diff(t_plot[peaks]) if len(peaks) > 1 else np.array([np.nan])
        mean_hr = 60 / np.mean(rr_intervals) if len(rr_intervals) > 0 else np.nan

        # Determine subplot indices
        if condition == "rest":
            raw_idx = 1
            filt_idx = 2
        else:  # exercise
            raw_idx = 3
            filt_idx = 4

        # --- Raw subplot ---
        plt.subplot(2, 2, raw_idx)
        plt.plot(t_plot, ecg_plot, color=color)
        plt.title(f"{condition.capitalize()} — Raw")
        plt.xlabel("Time [s]")
        plt.ylabel("Amplitude [mV]")
        plt.grid(True, linestyle=":", alpha=0.5)

        # --- Filtered subplot with R peaks ---
        plt.subplot(2, 2, filt_idx)
        plt.plot(t_plot, ecg_filt_plot, color=color)
        plt.plot(t_plot[peaks], ecg_filt_plot[peaks], 'ro', label='R peaks')
        plt.title(f"{condition.capitalize()} — Filtered\nMean HR: {mean_hr:.1f} BPM")
        plt.xlabel("Time [s]")
        plt.ylabel("Amplitude [mV]")
        plt.grid(True, linestyle=":", alpha=0.5)
        plt.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    save_path_ecg = os.path.join(OUTPUT_DIR, f"{subj}_ECG_Rpeaks.png")
    plt.savefig(save_path_ecg, dpi=300)
    plt.close()
    print(f"Saved ECG figure with R peaks: {save_path_ecg}")
