# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 15:19:19 2025

@author: sydne
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

# ---------- SETTINGS ----------
DATA_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Data"
OUTPUT_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Results/Raw V Filtered Data"
FS = 100  # Hz
PLOT_DURATION = 30 # seconds
OUTPUT_DIR = os.path.join(OUTPUT_DIR, f"output_{PLOT_DURATION}s")
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLORS = {
    "resting": "tab:blue",
    "rest": "tab:blue",
    "exercise": "tab:red",
}

# Ensure output folder exists
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
    if match:
        subj = match.group(1).upper()
    else:
        subj = "UNKNOWN"

    if subj not in subjects:
        subjects[subj] = []
    subjects[subj].append(file)

# ---------- PROCESS EACH SUBJECT ----------
for subj, subj_files in subjects.items():

    # ----- ECG FIGURE -----
    plt.figure(figsize=(12, 8))
    plt.suptitle(f"Subject {subj} — ECG Signals", fontsize=14, fontweight='bold')

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

        # Determine subplot index
        if condition == "rest":
            raw_idx = 1
            filt_idx = 2
        else:  # exercise
            raw_idx = 3
            filt_idx = 4

        # Raw
        plt.subplot(2, 2, raw_idx)
        plt.plot(t_plot, ecg_plot, color=color)
        plt.title(f"{condition.capitalize()} — Raw")
        plt.xlabel("Time [s]")
        plt.ylabel("Amplitude [mV]")
        plt.grid(True, linestyle=":", alpha=0.5)

        # Filtered
        plt.subplot(2, 2, filt_idx)
        plt.plot(t_plot, ecg_filt_plot, color=color)
        plt.title(f"{condition.capitalize()} — Filtered")
        plt.xlabel("Time [s]")
        plt.ylabel("Amplitude [mV]")
        plt.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    save_path_ecg = os.path.join(OUTPUT_DIR, f"{subj}_ECG.png")
    plt.savefig(save_path_ecg, dpi=300)
    plt.close()
    print(f"Saved ECG figure: {save_path_ecg}")

    # ----- RESPIRATION FIGURE -----
    plt.figure(figsize=(12, 8))
    plt.suptitle(f"Subject {subj} — Respiration Signals", fontsize=14, fontweight='bold')

    for file in subj_files:
        file_path = os.path.join(DATA_DIR, file)
        data = np.loadtxt(file_path, comments="#")
        resp = data[:, 6] * 1000  # convert to mV
        t = np.arange(len(resp)) / FS

        condition = "exercise" if "exercise" in file.lower() else "rest"
        color = COLORS.get(condition, "gray")
        resp_filt = bandpass_filter(resp, 0.05, 0.5, FS)

        # Limit to first 30 seconds
        num_samples = min(len(resp), FS * PLOT_DURATION)
        t_plot = t[:num_samples]
        resp_plot = resp[:num_samples]
        resp_filt_plot = resp_filt[:num_samples]

        # Determine subplot index
        if condition == "rest":
            raw_idx = 1
            filt_idx = 2
        else:
            raw_idx = 3
            filt_idx = 4

        # Raw
        plt.subplot(2, 2, raw_idx)
        plt.plot(t_plot, resp_plot, color=color)
        plt.title(f"{condition.capitalize()} — Raw")
        plt.xlabel("Time [s]")
        plt.ylabel("Amplitude [mV]")
        plt.grid(True, linestyle=":", alpha=0.5)

        # Filtered
        plt.subplot(2, 2, filt_idx)
        plt.plot(t_plot, resp_filt_plot, color=color)
        plt.title(f"{condition.capitalize()} — Filtered")
        plt.xlabel("Time [s]")
        plt.ylabel("Amplitude [mV]")
        plt.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    save_path_resp = os.path.join(OUTPUT_DIR, f"{subj}_Respiration.png")
    plt.savefig(save_path_resp, dpi=300)
    plt.close()
    print(f"Saved Respiration figure: {save_path_resp}")
