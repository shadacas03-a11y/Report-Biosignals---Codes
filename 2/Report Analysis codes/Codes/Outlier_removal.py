import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, find_peaks

# ---------- SETTINGS ----------
DATA_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Data"
OUTPUT_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Results/Outlier Removal"
PLOT_DURATION = 30  # seconds
FS = 100  # Hz

OUTPUT_DIR = os.path.join(OUTPUT_DIR, f"output_{PLOT_DURATION}s")
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

summary_ecg = {}
summary_resp = {}

for subj, subj_files in subjects.items():
    summary_ecg[subj] = {}
    summary_resp[subj] = {}

    # ---------- ECG FIGURE ----------
    fig_ecg, axes_ecg = plt.subplots(2, 2, figsize=(12, 8))
    fig_ecg.suptitle(f"Subject {subj} — ECG Signals with R Peaks", fontsize=14, fontweight='bold')

    # ---------- RESP FIGURE ----------
    fig_resp, axes_resp = plt.subplots(2, 2, figsize=(12, 8))
    fig_resp.suptitle(f"Subject {subj} — Respiration Signals with Peaks", fontsize=14, fontweight='bold')

    for file in subj_files:
        file_path = os.path.join(DATA_DIR, file)
        data = np.loadtxt(file_path, comments="#")
        condition = "exercise" if "exercise" in file.lower() else "rest"
        color = COLORS.get(condition, "gray")
        num_samples = min(len(data), FS * PLOT_DURATION)
        t_plot = np.arange(num_samples) / FS

        # ---------- ECG ----------
        ecg = data[:num_samples, 7] * 1000
        ecg_filt = bandpass_filter(ecg, 0.5, 40, FS)
        peaks_ecg = find_peaks(ecg_filt, distance=FS*0.5)[0]
        rr_intervals = np.diff(t_plot[peaks_ecg]) if len(peaks_ecg) > 1 else np.array([np.nan])
        rr_clean = remove_outliers(rr_intervals)
        mean_hr = 60 / np.mean(rr_clean) if len(rr_clean) > 0 else np.nan

        summary_ecg[subj][condition] = {
            "R_peaks": peaks_ecg,
            "RR_intervals": rr_intervals,
            "RR_intervals_clean": rr_clean,
            "mean_hr": mean_hr
        }

        # Determine subplot positions
        raw_idx = 0 if condition == "rest" else 1
        filt_idx = 2 if condition == "rest" else 3

        # Map subplot indices to axes
        axes_positions = {0: (0,0), 1: (1,0), 2: (0,1), 3: (1,1)}

        # Raw ECG
        ax = axes_ecg[axes_positions[raw_idx]]
        ax.plot(t_plot, ecg, color=color)
        ax.set_title(f"{condition.capitalize()} — Raw")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Amplitude [mV]")
        ax.grid(True, linestyle=":", alpha=0.5)

        # Filtered ECG
        ax = axes_ecg[axes_positions[filt_idx]]
        ax.plot(t_plot, ecg_filt, color=color)
        ax.plot(t_plot[peaks_ecg], ecg_filt[peaks_ecg], 'ro', label='R peaks')
        ax.set_title(f"{condition.capitalize()} — Filtered\nMean HR: {mean_hr:.1f} BPM")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Amplitude [mV]")
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.legend()

        # ---------- RESPIRATION ----------
        resp = data[:num_samples, 6] * 1000
        resp_filt = bandpass_filter(resp, 0.05, 0.5, FS)
        peaks_resp = find_peaks(resp_filt, distance=FS*1.0)[0]
        breath_intervals = np.diff(t_plot[peaks_resp]) if len(peaks_resp) > 1 else np.array([np.nan])
        breath_clean = remove_outliers(breath_intervals)
        mean_breath_rate = 60 / np.mean(breath_clean) if len(breath_clean) > 0 else np.nan

        summary_resp[subj][condition] = {
            "peaks": peaks_resp,
            "breath_intervals": breath_intervals,
            "breath_intervals_clean": breath_clean,
            "mean_breath_rate": mean_breath_rate
        }

        # Raw Respiration
        ax = axes_resp[axes_positions[raw_idx]]
        ax.plot(t_plot, resp, color=color)
        ax.set_title(f"{condition.capitalize()} — Raw")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Amplitude [mV]")
        ax.grid(True, linestyle=":", alpha=0.5)

        # Filtered Respiration
        ax = axes_resp[axes_positions[filt_idx]]
        ax.plot(t_plot, resp_filt, color=color)
        ax.plot(t_plot[peaks_resp], resp_filt[peaks_resp], 'ro', label='Peaks')
        ax.set_title(f"{condition.capitalize()} — Filtered\nMean BR: {mean_breath_rate:.1f} BPM")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Amplitude [mV]")
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.legend()

    # Save figures
    fig_ecg.tight_layout(rect=[0,0,1,0.95])
    fig_ecg.savefig(os.path.join(OUTPUT_DIR, f"{subj}_ECG_Rpeaks.png"), dpi=300)
    plt.close(fig_ecg)

    fig_resp.tight_layout(rect=[0,0,1,0.95])
    fig_resp.savefig(os.path.join(OUTPUT_DIR, f"{subj}_Respiration_Peaks.png"), dpi=300)
    plt.close(fig_resp)

    # ---------- Boxplots of cleaned intervals ----------
    conditions = list(summary_ecg[subj].keys())

    # ECG RR intervals
    rr_clean_data = [summary_ecg[subj][cond]["RR_intervals_clean"] for cond in conditions]
    plt.figure(figsize=(6,4))
    plt.boxplot(rr_clean_data, labels=conditions)
    plt.ylabel("RR Interval (s)")
    plt.title(f"{subj} — RR Interval Comparison (Outliers Removed)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{subj}_RR_Interval_Comparison_Cleaned.png"), dpi=300)
    plt.close()

    # Respiration intervals
    breath_clean_data = [summary_resp[subj][cond]["breath_intervals_clean"] for cond in conditions]
    plt.figure(figsize=(6,4))
    plt.boxplot(breath_clean_data, labels=conditions)
    plt.ylabel("Breath Interval (s)")
    plt.title(f"{subj} — Breath Interval Comparison (Outliers Removed)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{subj}_Breath_Interval_Comparison_Cleaned.png"), dpi=300)
    plt.close()
