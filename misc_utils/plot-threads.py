import matplotlib.pyplot as plt
import numpy as np

# Input data
threads = [1, 2, 4, 8, 16]
ram_usage_gb = [, , , , ]
execution_time_sec = [, , , , ]

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

plt.savefig('../results/threads-stats-vg.png', format='png')
plt.savefig('../results/threads-stats-vg.eps', format='eps')