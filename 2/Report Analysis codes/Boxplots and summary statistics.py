# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 16:42:25 2025

@author: sydne
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks

# ---------- SETTINGS ----------
DATA_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Data"
OUTPUT_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Results/Boxplot Results"
PLOT_DURATION = 30  # seconds
FS = 100  # Hz

OUTPUT_DIR = os.path.join(OUTPUT_DIR, f"output_{PLOT_DURATION}s")
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLORS = {"rest": "tab:blue", "exercise": "tab:red"}

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
    q1, q3 = np.percentile(data, [25, 75])
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    return data[(data >= lower) & (data <= upper)]

# ---------- GROUP FILES BY SUBJECT ----------
files = [f for f in os.listdir(DATA_DIR) if f.endswith(".txt")]
subjects = {}
for file in files:
    match = re.match(r"([AB])_", file, re.IGNORECASE)
    subj = match.group(1).upper() if match else "UNKNOWN"
    subjects.setdefault(subj, []).append(file)

# ---------- MASTER STATS LIST ----------
all_stats = []

# ---------- MAIN LOOP ----------
for subj, subj_files in subjects.items():
    print(f"\n🔹 Processing Subject {subj}")
    summary_ecg = {}
    summary_resp = {}

    for file in subj_files:
        file_path = os.path.join(DATA_DIR, file)
        data = np.loadtxt(file_path, comments="#")

        condition = "exercise" if "exercise" in file.lower() else "rest"
        num_samples = min(len(data), FS * PLOT_DURATION)
        t_plot = np.arange(num_samples) / FS

        # ---------- ECG ----------
        ecg = data[:num_samples, 7] * 1000  # mV
        ecg_filt = bandpass_filter(ecg, 0.5, 40, FS)
        peaks_ecg, _ = find_peaks(ecg_filt, distance=FS*0.5)
        rr_intervals = np.diff(t_plot[peaks_ecg]) if len(peaks_ecg) > 1 else np.array([])
        rr_clean = remove_outliers(rr_intervals)
        mean_hr = 60 / np.mean(rr_clean) if len(rr_clean) > 0 else np.nan

        summary_ecg[condition] = {
            "RR_intervals": rr_intervals,
            "RR_intervals_clean": rr_clean,
            "mean_hr": mean_hr
        }

        # ---------- RESPIRATION ----------
        resp = data[:num_samples, 6] * 1000  # mV
        resp_filt = bandpass_filter(resp, 0.05, 0.5, FS)
        peaks_resp, _ = find_peaks(resp_filt, distance=FS*1.0)
        breath_intervals = np.diff(t_plot[peaks_resp]) if len(peaks_resp) > 1 else np.array([])
        breath_clean = remove_outliers(breath_intervals)
        mean_br = 60 / np.mean(breath_clean) if len(breath_clean) > 0 else np.nan

        summary_resp[condition] = {
            "breath_intervals": breath_intervals,
            "breath_intervals_clean": breath_clean,
            "mean_breath_rate": mean_br
        }

    # ---------- BOXPLOTS ----------
    for signal_type, summary in [("RR", summary_ecg), ("Resp", summary_resp)]:
        raw_data, clean_data = [], []

        for cond, d in summary.items():
            if signal_type == "RR":
                raw = d["RR_intervals"]
                clean = d["RR_intervals_clean"]
                ylabel = "RR Interval (s)"
            else:
                raw = d["breath_intervals"]
                clean = d["breath_intervals_clean"]
                ylabel = "Respiration Interval (s)"

            if len(raw) > 0:
                raw_data.append(pd.DataFrame({"Interval": raw, "Condition": cond, "Type": "Before"}))
            if len(clean) > 0:
                clean_data.append(pd.DataFrame({"Interval": clean, "Condition": cond, "Type": "After"}))

        if len(raw_data) == 0 and len(clean_data) == 0:
            continue

        df_combined = pd.concat(raw_data + clean_data)
        plt.figure(figsize=(8,6))
        sns.boxplot(data=df_combined, x="Condition", y="Interval", hue="Type", palette="Set2")
        plt.title(f"{subj} — {signal_type}-intervals Before and After Outlier Removal")
        plt.ylabel(ylabel)
        plt.xlabel("Condition")
        plt.legend(title="Outlier Removal Stage")
        plt.tight_layout()
        save_path = os.path.join(OUTPUT_DIR, f"{subj}_{signal_type}_Boxplot_BeforeAfter.png")
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"✅ Saved boxplot: {save_path}")

    # ---------- COLLECT SUMMARY STATISTICS ----------
    def compute_stats(data, subject, cond, signal_type, stage):
        if len(data) == 0:
            return None
        return {
            "Subject": subject,
            "Condition": cond,
            "Signal": signal_type,
            "Stage": stage,
            "Mean": np.mean(data),
            "Median": np.median(data),
            "Std": np.std(data),
            "Min": np.min(data),
            "Max": np.max(data),
            "Count": len(data)
        }

    for cond, d in summary_ecg.items():
        all_stats.append(compute_stats(d["RR_intervals"], subj, cond, "RR", "Before"))
        all_stats.append(compute_stats(d["RR_intervals_clean"], subj, cond, "RR", "After"))

    for cond, d in summary_resp.items():
        all_stats.append(compute_stats(d["breath_intervals"], subj, cond, "Resp", "Before"))
        all_stats.append(compute_stats(d["breath_intervals_clean"], subj, cond, "Resp", "After"))

# ---------- SAVE MASTER SUMMARY CSV ----------
df_stats = pd.DataFrame([s for s in all_stats if s is not None])
csv_path = os.path.join(OUTPUT_DIR, "All_Subjects_Summary_Statistics.csv")
df_stats.to_csv(csv_path, index=False)
print(f"\n📊 Saved master summary statistics CSV: {csv_path}")
