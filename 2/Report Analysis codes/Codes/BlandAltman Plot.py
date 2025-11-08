# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 16:49:11 2025

@author: sydne
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ========================
# CONFIGURATION
# ========================
DATA_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Results/Bland Altman Plot Data"
CSV_FILE = os.path.join(DATA_DIR, "combined_hr_summary.csv")
OUTPUT_DIR = os.path.join(DATA_DIR, "bland_altman_plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ========================
# LOAD COMBINED HR DATA
# ========================
df = pd.read_csv(CSV_FILE)

# ========================
# FUNCTION TO CREATE BLAND-ALTMAN PLOT
# ========================
def bland_altman_plot(hr_ecg, hr_pulse, subject, condition, save_path):
    mean_hr = (hr_ecg + hr_pulse) / 2
    diff_hr = hr_pulse - hr_ecg

    mean_diff = np.nanmean(diff_hr)
    sd_diff = np.nanstd(diff_hr)
    loa_upper = mean_diff + 1.96 * sd_diff
    loa_lower = mean_diff - 1.96 * sd_diff

    plt.figure(figsize=(8,6))
    plt.scatter(mean_hr, diff_hr, alpha=0.7, edgecolor='k', label=f'{subject} {condition}')
    plt.axhline(mean_diff, color='blue', linestyle='--', label=f'Mean Bias = {mean_diff:.2f}')
    plt.axhline(loa_upper, color='red', linestyle='--', label=f'+1.96 SD = {loa_upper:.2f}')
    plt.axhline(loa_lower, color='red', linestyle='--', label=f'-1.96 SD = {loa_lower:.2f}')

    plt.title(f'Bland–Altman Plot — {subject} {condition.capitalize()}')
    plt.xlabel('Mean HR (BPM)')
    plt.ylabel('Difference (Pulse − ECG) (BPM)')
    plt.legend(loc='upper right')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"✅ Saved plot to {save_path}")

# ========================
# CREATE PLOTS FOR EACH CASE
# ========================
for subject in df['Subject'].unique():
    for condition in df['Condition'].unique():
        subset = df[(df['Subject'] == subject) & (df['Condition'] == condition)]
        if subset.empty:
            continue
        save_path = os.path.join(OUTPUT_DIR, f'BlandAltman_{subject}_{condition}.png')
        bland_altman_plot(
            subset['Avg_HR_ECG_BPM'].values,
            subset['Avg_HR_Pulse_BPM'].values,
            subject,
            condition,
            save_path
        )

print("/nAll Bland–Altman plots saved.")
