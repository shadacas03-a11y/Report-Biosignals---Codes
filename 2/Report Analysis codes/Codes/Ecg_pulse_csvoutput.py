import os
import re
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks

# ========================
# CONFIGURATION
# ========================
DATA_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Data"
FS = 100                               # sampling frequency (Hz)
time_windows = [(0,30), (30,60)]       # list of (start_sec, end_sec)
OUTPUT_DIR = "C:/Users/sydne/OneDrive/Desktop/MScTE Env M1/Biosignals/2. Consequences of exercising/Results/Bland Altman Plot Data"
os.makedirs(OUTPUT_DIR, exist_ok=True)
SUMMARY_FILE = os.path.join(OUTPUT_DIR, "combined_hr_summary.csv")

FS = 100  # sampling rate (Hz)
WINDOW_DURATION = 30  # seconds

# Filter settings
ECG_LOW, ECG_HIGH = 0.5, 40   # ECG bandpass range (Hz)
PULSE_LOW, PULSE_HIGH = 0.5, 5  # Pulse bandpass range (Hz)

# ========================
# HELPER FUNCTIONS
# ========================
def butter_bandpass(lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low, high = lowcut / nyq, highcut / nyq
    b, a = butter(order, [low, high], btype="band")
    return b, a

def bandpass_filter(data, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order)
    return filtfilt(b, a, data)

def compute_hr_from_signal(signal, fs, is_ecg=True):
    """Compute average HR (BPM) for a given signal window."""
    if is_ecg:
        # ECG: sharper peaks, more frequent → stricter distance
        peaks, _ = find_peaks(signal, distance=fs * 0.4, height=np.percentile(signal, 75))
    else:
        # Pulse: smoother waveform → relaxed threshold
        peaks, _ = find_peaks(signal, distance=fs * 0.4, height=np.percentile(signal, 70))

    if len(peaks) < 2:
        return np.nan

    rr_intervals = np.diff(peaks) / fs
    bpm = 60 / np.mean(rr_intervals)
    return bpm

# ========================
# MAIN SCRIPT
# ========================
summary_records = []
files = [f for f in os.listdir(DATA_DIR) if f.endswith(".txt")]

for file in files:
    filepath = os.path.join(DATA_DIR, file)
    print(f"Processing {file}...")

    # Identify subject and condition from filename
    subj_match = re.search(r"([AB])_", file)
    cond_match = re.search(r"(rest|exercise)", file.lower())
    subject = subj_match.group(1) if subj_match else "Unknown"
    condition = cond_match.group(1) if cond_match else "unspecified"

    # Load raw file
    df = pd.read_csv(filepath, comment="#", sep="\t", header=None)
    df = df.dropna(axis=1, how="all")

    # Extract relevant channels
    pulse_raw = df.iloc[:, 5].values  # A2 = pulse
    ecg_raw = df.iloc[:, 7].values    # A4 = ECG

    # Convert to millivolts (assuming 10-bit ADC, 3.3 V reference)
    pulse_mV = (pulse_raw / 1023.0) * 3.3 * 1000
    ecg_mV = (ecg_raw / 1023.0) * 3.3 * 1000

    # Filter both signals
    pulse_filt = bandpass_filter(pulse_mV, PULSE_LOW, PULSE_HIGH, FS)
    ecg_filt = bandpass_filter(ecg_mV, ECG_LOW, ECG_HIGH, FS)

    # Segment into 30 s windows
    total_duration = len(ecg_filt) / FS
    num_windows = int(total_duration // WINDOW_DURATION)

    for i in range(num_windows):
        start = int(i * WINDOW_DURATION * FS)
        end = int((i + 1) * WINDOW_DURATION * FS)

        pulse_segment = pulse_filt[start:end]
        ecg_segment = ecg_filt[start:end]

        hr_pulse = compute_hr_from_signal(pulse_segment, FS, is_ecg=False)
        hr_ecg = compute_hr_from_signal(ecg_segment, FS, is_ecg=True)

        summary_records.append({
            "Subject": subject,
            "Condition": condition,
            "Window_Number": i + 1,
            "Start_Time_s": start / FS,
            "End_Time_s": end / FS,
            "Avg_HR_Pulse_BPM": round(hr_pulse, 2),
            "Avg_HR_ECG_BPM": round(hr_ecg, 2),
        })

# ========================
# SAVE RESULTS
# ========================
summary_df = pd.DataFrame(summary_records)
summary_df.to_csv(SUMMARY_FILE, index=False)

print(f"\n✅ Combined HR summary saved to: {SUMMARY_FILE}")
print(summary_df.head(10))
