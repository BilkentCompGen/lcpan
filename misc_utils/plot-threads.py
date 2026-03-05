import matplotlib.pyplot as plt
import numpy as np

# Input data
threads = [1, 2, 4, 8, 16]
lcpan_ram_usage_gb = [6.99, 7.45, 7.54, 8.20, 8.31]
lcpan_execution_time_sec = [340.078, 219.522, 168.260, 155.100, 162.424]
vg_ram_usage_gb = [114.43, 114.42, 114.46, 114.40, 114.37]
vg_execution_time_sec = [5059.42, 4894.38, 4838.08, 4969.56, 4925.21]
vg_gnu_ram_usage_gb = [114.53, 114.53, 114.40, 114.55, 114.52]
vg_gnu_execution_time_sec = [4954.81, 3821.73, 3309.53, 3008.64, 2929.77]

x = np.arange(len(threads))
width = 0.30

# ---- Create broken-axis layout ----
fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, sharex=True, figsize=(9, 7),
    gridspec_kw={'height_ratios': [2, 2], 'hspace': 0.05}
)

# Twin axes for execution time
ax_top_time = ax_top.twinx()
ax_bot_time = ax_bot.twinx()

# ---- Bottom (LOW VALUES) ----
ax_bot.bar(x - width, lcpan_ram_usage_gb, width,
           color='#aec7e8', alpha=0.7, label='LCPan RAM')
ax_bot.bar(x, vg_ram_usage_gb, width,
           color='#ffbb78', alpha=0.7, label='VG RAM')
ax_bot.bar(x + width, vg_gnu_ram_usage_gb, width,
           color='#98df8a', alpha=0.7, label='VG-GNU RAM')

ax_bot_time.plot(x - width, lcpan_execution_time_sec,
                 color='#1f77b4', marker='o', linestyle='--', label='LCPan Time')
ax_bot_time.plot(x, vg_execution_time_sec,
                 color='#ff7f0e', marker='o', linestyle='--', label='VG Time')
ax_bot_time.plot(x + width, vg_gnu_execution_time_sec,
                 color='#2ca02c', marker='o', linestyle='--', label='VG-GNU Time')

ax_bot.set_ylim(0, 9)
ax_bot_time.set_ylim(0, 500)

# ---- Top (HIGH VALUES) ----
ax_top.bar(x - width, lcpan_ram_usage_gb, width,
           color='#aec7e8', alpha=0.7)
ax_top.bar(x, vg_ram_usage_gb, width,
           color='#ffbb78', alpha=0.7)
ax_top.bar(x + width, vg_gnu_ram_usage_gb, width,
           color='#98df8a', alpha=0.7)

ax_top_time.plot(x - width, lcpan_execution_time_sec,
                 color='#1f77b4', marker='o', linestyle='--')
ax_top_time.plot(x, vg_execution_time_sec,
                 color='#ff7f0e', marker='o', linestyle='--')
ax_top_time.plot(x + width, vg_gnu_execution_time_sec,
                 color='#2ca02c', marker='o', linestyle='--')

ax_top.set_ylim(110, 116)
ax_top_time.set_ylim(2700, 5600)

# ---- Labels ----
ax_bot.set_xlabel('# Threads')
ax_bot.set_ylabel('')
ax_top.set_ylabel('')
ax_bot_time.set_ylabel('')
ax_top_time.set_ylabel('')

fig.text(0.03, 0.5, 'RAM Usage (GB)',
         va='center', rotation='vertical')

fig.text(0.97, 0.5, 'Execution Time (s)',
         va='center', rotation='vertical', ha='right')

ax_bot.set_xticks(x)
ax_bot.set_xticklabels(threads)

# ---- Break marks ----
# ---- Subtle break marks ----
d = 0.008
kwargs = dict(color='k', clip_on=False, linewidth=0.8)

ax_top.plot((-d, +d), (-d, +d), transform=ax_top.transAxes, **kwargs)
ax_top.plot((1-d, 1+d), (-d, +d), transform=ax_top.transAxes, **kwargs)

ax_bot.plot((-d, +d), (1-d, 1+d), transform=ax_bot.transAxes, **kwargs)
ax_bot.plot((1-d, 1+d), (1-d, 1+d), transform=ax_bot.transAxes, **kwargs)

# ---- Legend ----
handles1, labels1 = ax_bot.get_legend_handles_labels()
handles2, labels2 = ax_bot_time.get_legend_handles_labels()

handles = []
labels = []

for h1, h2, l1, l2 in zip(handles1, handles2, labels1, labels2):
    handles.extend([h1, h2])
    labels.extend([l1, l2])

ax_top.legend(handles, labels, ncol=3, loc='upper left')

plt.subplots_adjust(
    left=0.09,
    right=0.89,
    top=0.96,
    bottom=0.10,
    hspace=0.05
)

plt.savefig('../results/threads-stats-vg.png', format='png')
plt.savefig('../results/threads-stats-vg.eps', format='eps')