import matplotlib.pyplot as plt
import numpy as np

# Input data
threads = [1, 2, 4, 8, 16]
ram_usage_kb = [7250784, 7790140, 7790696, 8809136, 8945664]  # KB
execution_time_str = ["12:24.05", "8:12.98", "5:57.42", "4:48.25", "5:02.95"]  # H:MM:SS.ss or MM:SS.ss

# Convert RAM usage to GB
ram_usage_gb = [x / (1024 ** 2) for x in ram_usage_kb]  # 1 GB = 1024^2 KB

# Convert execution time to seconds
def time_to_seconds(time_str):
    parts = time_str.split(":")
    if len(parts) == 3:  # H:MM:SS.ss
        hours, minutes, seconds = int(parts[0]), int(parts[1]), float(parts[2])
    elif len(parts) == 2:  # MM:SS.ss
        hours, minutes, seconds = 0, int(parts[0]), float(parts[1])
    else:
        raise ValueError(f"Invalid time format: {time_str}")
    return hours * 3600 + minutes * 60 + seconds

execution_time_sec = [time_to_seconds(t) for t in execution_time_str]

fig, ax1 = plt.subplots()

# Ensure even spacing on X-axis
ax1.set_xticks(range(len(threads)))
ax1.set_xticklabels(threads)

# Bar plot for RAM usage
ax1.bar(range(len(threads)), ram_usage_gb, color='royalblue', alpha=0.7, label='RAM Usage (GB)')
ax1.set_xlabel('Threads')
ax1.set_ylabel('RAM Usage (GB)', color='royalblue')
ax1.tick_params(axis='y', labelcolor='royalblue')
ax1.set_ylim(0, max(ram_usage_gb) * 1.1)  # Start from 0 and add margin

# Second Y-axis for Execution Time
ax2 = ax1.twinx()
ax2.plot(range(len(threads)), execution_time_sec, color='darkorange', marker='o', linestyle='-', linewidth=2, label='Execution Time (s)')
ax2.set_ylabel('Execution Time (s)', color='darkorange')
ax2.tick_params(axis='y', labelcolor='darkorange')
ax2.set_ylim(0, max(execution_time_sec) * 1.1)  # Start from 0 and add margin

# Grid and Layout
ax1.grid(True, linestyle='--', alpha=0.6)
fig.subplots_adjust(top=0.90, bottom=0.10)  # Increase margins

# Save Plot
# plt.title('RAM Usage and Execution Time vs Threads')
plt.savefig('threads-stats-lcpan.eps', format='eps')
