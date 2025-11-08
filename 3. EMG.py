#!/usr/bin/env python
# coding: utf-8

# In[8]:


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

# -----------------------------
# Parameters
# -----------------------------
folder_path = r"C:\Users\rodri\EMG data processing"  
fs = 1000
lowcut = 20
highcut = 450

files_info = {
    "studentA_flexion.txt": ("Student A", "Flexion-Extension"),
    "studentB_flexion.txt": ("Student B", "Flexion-Extension"),
    "studentA-pushups.txt": ("Student A", "Push-Ups"),
    "studentB-pushups.txt": ("Student B", "Push-Ups")
}

# Color scheme for student+activity
color_scheme = {
    ("Student A", "Flexion-Extension"): ("blue", "cyan"),
    ("Student B", "Flexion-Extension"): ("red", "orange"),
    ("Student A", "Push-Ups"): ("green", "lime"),
    ("Student B", "Push-Ups"): ("purple", "magenta")
}

# -----------------------------
# Functions
# -----------------------------
def read_emg_file(file_path):
    """Read EMG file from OpenSignals and return A3 and A4."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    # Find header end
    header_end_idx = 0
    for i, line in enumerate(lines):
        if "# EndOfHeader" in line:
            header_end_idx = i
            break
    df = pd.read_csv(file_path, sep='\t', skiprows=header_end_idx+1, header=None,
                     usecols=range(7),  # keep only first 7 columns
                     names=["nSeq","I1","I2","O1","O2","A3","A4"])
    # Print first 5 rows to check
    print(f"\nFirst 5 rows of {os.path.basename(file_path)}:")
    print(df.head())
    return df[["A3","A4"]]

def bandpass_filter(data, fs, lowcut, highcut, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    if len(data) <= 3 * max(len(a), len(b)):
        print("Warning: data too short for filtfilt, skipping filter")
        return data
    return filtfilt(b, a, data)

def plot_emg_signal(emg_df, student, activity, save_path):
    t = np.arange(len(emg_df)) / fs
    color_a3, color_a4 = color_scheme[(student, activity)]

    plt.figure(figsize=(12,5))
    plt.plot(t, emg_df["A3"], label="A3 (Biceps)", color=color_a3)
    plt.plot(t, emg_df["A4"], label="A4 (Triceps)", color=color_a4)
    plt.title(f"EMG Signal - {student} - {activity}")
    plt.xlabel("Time (s)")
    plt.ylabel("EMG (V)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

# -----------------------------
# Main Processing
# -----------------------------
for file_name, (student, activity) in files_info.items():
    file_path = os.path.join(folder_path, file_name)

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        continue

    # Read EMG and check first 5 rows
    emg_data = read_emg_file(file_path)

    # Filter
    emg_filtered = emg_data.apply(lambda x: bandpass_filter(x.values, fs, lowcut, highcut))

    # Plot and save
    save_file = os.path.join(folder_path, f"{file_name.replace('.txt','')}_filtered.png")
    plot_emg_signal(emg_filtered, student, activity, save_file)

    print(f"Processed and saved: {save_file}")


# In[9]:


# -----------------------------
# Function to full-wave rectify
# -----------------------------
def rectify_emg(emg_df):
    """Full-wave rectify EMG by taking the absolute value."""
    return emg_df.abs()

# -----------------------------
# Main Processing (Step 2)
# -----------------------------
for file_name, (student, activity) in files_info.items():
    file_path = os.path.join(folder_path, file_name)

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        continue

    # Read EMG
    emg_data = read_emg_file(file_path)

    # Bandpass filter
    emg_filtered = emg_data.apply(lambda x: bandpass_filter(x.values, fs, lowcut, highcut))

    # Full-wave rectify
    emg_rectified = rectify_emg(emg_filtered)

    # Plot and save
    save_file = os.path.join(folder_path, f"{file_name.replace('.txt','')}_rectified.png")
    plot_emg_signal(emg_rectified, student, activity, save_file)

    print(f"Processed and saved (rectified): {save_file}")


# In[12]:


import matplotlib.pyplot as plt

# -----------------------------
# Step: Compute and store MVC
# -----------------------------
mvc_values = {}

for student in ["Student A", "Student B"]:
    file_name = f"{student.lower().replace(' ','')}_flexion.txt"
    file_path = os.path.join(folder_path, file_name)

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        continue

    # Read, filter, rectify
    emg_data = read_emg_file(file_path)
    emg_filtered = emg_data.apply(lambda x: bandpass_filter(x.values, fs, lowcut, highcut))
    emg_rectified = rectify_emg(emg_filtered)

    # Compute MVC for each muscle
    mvc_values[student] = {
        "A3": emg_rectified["A3"].max(),
        "A4": emg_rectified["A4"].max()
    }
    print(f"{student} MVC values: {mvc_values[student]}")

# -----------------------------
# Step: Plot MVC bar chart
# -----------------------------
students = list(mvc_values.keys())
muscles = ["A3 (Biceps)", "A4 (Triceps)"]

# Extract values for plotting
mvc_biceps = [mvc_values[s]["A3"] for s in students]
mvc_triceps = [mvc_values[s]["A4"] for s in students]

x = np.arange(len(students))  # positions
width = 0.35

plt.figure(figsize=(8,5))
plt.bar(x - width/2, mvc_biceps, width, label='Biceps (A3)', color='blue')
plt.bar(x + width/2, mvc_triceps, width, label='Triceps (A4)', color='red')
plt.xticks(x, students)
plt.ylabel("MVC (V)")
plt.title("Maximum Voluntary Contraction (MVC) per Student")
plt.legend()
plt.tight_layout()

# Save figure
mvc_plot_path = os.path.join(folder_path, "MVC_values.png")
plt.savefig(mvc_plot_path)
plt.show()
print(f"MVC bar chart saved as: {mvc_plot_path}")


# In[13]:


from scipy.signal import butter, filtfilt

# -----------------------------
# Low-pass filter for envelope
# -----------------------------
def lowpass_filter(data, fs, cutoff=5, order=4):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low')
    return filtfilt(b, a, data)

# -----------------------------
# Step 4: Plot flexion-extension for each student
# -----------------------------
for student in ["Student A", "Student B"]:
    file_name = f"{student.lower().replace(' ','')}_flexion.txt"
    file_path = os.path.join(folder_path, file_name)

    if not os.path.exists(file_path):
        continue

    # Read, filter, rectify, normalize
    emg_data = read_emg_file(file_path)
    emg_filtered = emg_data.apply(lambda x: bandpass_filter(x.values, fs, lowcut, highcut))
    emg_rectified = rectify_emg(emg_filtered)
    emg_normalized = normalize_by_mvc(emg_rectified, mvc_values[student])

    # Smooth to get envelope
    emg_envelope = emg_normalized.apply(lambda x: lowpass_filter(x.values, fs, cutoff=5))

    # Plot
    t = np.arange(len(emg_envelope)) / fs
    color_a3, color_a4 = color_scheme[(student, "Flexion-Extension")]

    plt.figure(figsize=(12,5))
    plt.plot(t, emg_envelope["A3"], label="Biceps (A3)", color=color_a3)
    plt.plot(t, emg_envelope["A4"], label="Triceps (A4)", color=color_a4)
    plt.title(f"Flexion-Extension EMG Envelope - {student}")
    plt.xlabel("Time (s)")
    plt.ylabel("Normalized EMG (0-1)")
    plt.legend()
    plt.tight_layout()

    save_file = os.path.join(folder_path, f"{student}_flexion_envelope.png")
    plt.savefig(save_file)
    plt.show()
    print(f"Saved envelope plot: {save_file}")


# In[15]:


# -----------------------------
# Step 5: Push-ups plot only
# -----------------------------
for student in ["Student A", "Student B"]:
    file_name = f"{student.lower().replace(' ','')}-pushups.txt"
    file_path = os.path.join(folder_path, file_name)

    if not os.path.exists(file_path):
        continue

    # Read EMG
    emg_data = read_emg_file(file_path)

    # Bandpass filter
    emg_filtered = emg_data.apply(lambda x: bandpass_filter(x.values, fs, lowcut, highcut))

    # Full-wave rectify
    emg_rectified = rectify_emg(emg_filtered)

    # Normalize by MVC
    emg_normalized = normalize_by_mvc(emg_rectified, mvc_values[student])

    # Smooth to get envelope
    emg_envelope = emg_normalized.apply(lambda x: lowpass_filter(x.values, fs, cutoff=5))

    # Plot EMG envelope
    t = np.arange(len(emg_envelope)) / fs
    color_a3, color_a4 = color_scheme[(student, "Push-Ups")]

    plt.figure(figsize=(12,5))
    plt.plot(t, emg_envelope["A3"], label="Biceps (A3)", color=color_a3)
    plt.plot(t, emg_envelope["A4"], label="Triceps (A4)", color=color_a4)
    plt.title(f"Push-Ups EMG Envelope - {student}")
    plt.xlabel("Time (s)")
    plt.ylabel("Normalized EMG (0-1)")
    plt.legend()
    plt.tight_layout()

    save_file = os.path.join(folder_path, f"{student}_pushups_envelope.png")
    plt.savefig(save_file)
    plt.show()
    print(f"Saved push-ups plot: {save_file}")


# In[16]:


# -----------------------------
# Step 7: Co-activation ratio
# -----------------------------
co_activation_results = {}

for student in ["Student A", "Student B"]:
    file_name = f"{student.lower().replace(' ','')}-pushups.txt"
    file_path = os.path.join(folder_path, file_name)

    if not os.path.exists(file_path):
        continue

    # Read, filter, rectify, normalize
    emg_data = read_emg_file(file_path)
    emg_filtered = emg_data.apply(lambda x: bandpass_filter(x.values, fs, lowcut, highcut))
    emg_rectified = rectify_emg(emg_filtered)
    emg_normalized = normalize_by_mvc(emg_rectified, mvc_values[student])

    # Smooth for envelope
    emg_envelope = emg_normalized.apply(lambda x: lowpass_filter(x.values, fs, cutoff=5))

    # Compute co-activation ratio
    EMG_biceps = emg_envelope["A3"].values
    EMG_triceps = emg_envelope["A4"].values
    Co = EMG_biceps / (EMG_biceps + EMG_triceps + 1e-8)  # avoid division by zero
    Co_avg = np.mean(Co)

    co_activation_results[student] = Co_avg
    print(f"{student} co-activation ratio (push-ups): {Co_avg:.3f}")


# In[17]:


import os
import numpy as np

# files and students (use exact filenames you have)
files_info = {
    "studentA_flexion.txt": ("Student A", "Flexion-Extension"),
    "studentB_flexion.txt": ("Student B", "Flexion-Extension"),
    "studentA-pushups.txt": ("Student A", "Push-Ups"),
    "studentB-pushups.txt": ("Student B", "Push-Ups")
}

# print existing MVC (from your pipeline) and push-up maxima
print("=== MVC and Push-up maxima diagnostics ===")
for student in ["Student A", "Student B"]:
    # print MVC you computed earlier (if stored in mvc_values)
    mvc = mvc_values.get(student, None)
    print(f"\n{student} stored MVC: {mvc}")

    # load push-up trial, process same as pipeline (filter->rectify->envelope)
    push_file = [fn for fn,info in files_info.items() if info[0]==student and "pushups" in fn.lower()]
    if not push_file:
        print("  push-up file not found in files_info")
        continue
    push_path = os.path.join(folder_path, push_file[0])
    if not os.path.exists(push_path):
        print(f"  push-up file not found on disk: {push_path}")
        continue

    emg_data = read_emg_file(push_path)
    emg_filtered = emg_data.apply(lambda x: bandpass_filter(x.values, fs, lowcut, highcut))
    emg_rectified = rectify_emg(emg_filtered)
    emg_envelope = emg_rectified.apply(lambda x: lowpass_filter(x.values, fs, cutoff=5))

    max_biceps_push = emg_rectified["A3"].max()
    max_triceps_push = emg_rectified["A4"].max()
    env_max_biceps = emg_envelope["A3"].max()
    env_max_triceps = emg_envelope["A4"].max()

    print(f"  Push-ups raw-rectified max: A3={max_biceps_push:.4f}, A4={max_triceps_push:.4f}")
    print(f"  Push-ups envelope max:     A3={env_max_biceps:.4f}, A4={env_max_triceps:.4f}")

    if mvc is not None:
        print(f"  Ratio (push max / MVC): A3={max_biceps_push/mvc['A3']:.3f}, A4={max_triceps_push/mvc['A4']:.3f}")
    else:
        print("  No stored MVC for this student.")


# In[ ]:




