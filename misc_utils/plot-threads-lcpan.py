import matplotlib.pyplot as plt
import numpy as np

# Input data
threads = [1, 2, 4, 8, 16]
# LCPan VG
ram_usage_kb = [0, 0, 0, 0, 0]  # KB
execution_time_str = ["", "", "", "", ""]  # H:MM:SS.ss or MM:SS.ss
merge_time_str = ["", "", "", "", ""]  # H:MM:SS.ss or MM:SS.ss
# LCPan VGX
# ram_usage_kb = [0, 0, 0, 0, 0]  # KB
# execution_time_str = ["", "", "", "", ""]  # H:MM:SS.ss or MM:SS.ss
# merge_time_str = ["", "", "", "", ""]  # H:MM:SS.ss or MM:SS.ss
# VG-toolkit
# ram_usage_kb = [0, 0, 0, 0, 0]  # KB
# execution_time_str = ["", "", "", "", ""]  # H:MM:SS.ss or MM:SS.ss
# merge_time_str = ["", "", "", "", ""]

# Convert RAM usage to GB
ram_usage_gb = [x / (1024 ** 2) for x in ram_usage_kb]  # 1 GB = 1024^2 KB

# Convert execution time to seconds
def time_to_seconds(exec_str, merge_str):
    exec_parts = exec_str.split(":")
    merge_parts = merge_str.split(":")
    if len(exec_parts) == 3:  # H:MM:SS.ss
        exec_hours, exec_minutes, exec_seconds = int(exec_parts[0]), int(exec_parts[1]), float(exec_parts[2])
    elif len(exec_parts) == 2:  # MM:SS.ss
        exec_hours, exec_minutes, exec_seconds = 0, int(exec_parts[0]), float(exec_parts[1])
    else:
        raise ValueError(f"Invalid time format: {exec_str}")
    if len(merge_parts) == 3:  # H:MM:SS.ss
        merge_hours, merge_minutes, merge_seconds = int(merge_parts[0]), int(merge_parts[1]), float(merge_parts[2])
    elif len(merge_parts) == 2:  # MM:SS.ss
        merge_hours, merge_minutes, merge_seconds = 0, int(merge_parts[0]), float(merge_parts[1])
    else:
        raise ValueError(f"Invalid time format: {merge_str}")
    return (exec_hours * 3600 + exec_minutes * 60 + exec_seconds) + (merge_hours * 3600 + merge_minutes * 60 + merge_seconds)

execution_time_sec = []
for i in range(len(execution_time_str)):
    execution_time_sec.append(time_to_seconds(execution_time_str[i], merge_time_str[i]))

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

plt.savefig('threads-stats-lcpan-vg.png', format='png')
# plt.savefig('threads-stats-lcpan-vgx.png', format='png')
# plt.savefig('threads-stats-vg.png', format='png')

# plt.savefig('threads-stats.eps', format='eps')