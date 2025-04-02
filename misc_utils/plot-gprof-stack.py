import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

data_t1 = """""".strip().split("\n")

data_t2 = """""".strip().split("\n")

data_t4 = """""".strip().split("\n")

data_t8 = """""".strip().split("\n")

data_t16 = """""".strip().split("\n")

all_data = [data_t1, data_t2, data_t4, data_t8, data_t16]

categories = ["1", "2", "4", "8", "16"]
data = {}
functions = []

for t in all_data:
    for line in t:
        func = line.split()[-1].strip()
        if func not in functions:
            functions.append(func)

# print(functions)

functions_to_be_printed = ['process_chrom', 'find_boundaries', 'variate', 'variate_snp', 'variate_sv', 'print_ref_seq', 'print_seq', 'print_link', 't_read_vcf', 'line_queue_push', 'line_queue_pop']

for t in all_data:
    for i, line in enumerate(t):
        values = line.split()
        function = values[-1].strip()
        if function not in data:
            data[function] = [0] * len(categories)
        data[function][all_data.index(t)] = float(values[0])

# Plotting
fig, ax = plt.subplots(figsize=(8, 6))

bottom = np.zeros(len(categories))
cmap = plt.colormaps.get_cmap('tab20')
colors = [cmap(i / len(functions_to_be_printed)) for i in range(len(functions_to_be_printed))] 

# Plot the stacked bars
for idx, function in enumerate(functions_to_be_printed):
    ax.bar(categories, data[function], label=function, bottom=bottom, color=colors[idx])
    bottom += np.array(data[function])

all_values = np.sum([np.array(data[func]) for func in data], axis=0)
selected_values = np.sum([np.array(data[func]) for func in functions_to_be_printed], axis=0)
other_values = all_values - selected_values

# Plot "Other" if it has nonzero values
if np.any(other_values > 0):
    ax.bar(categories, other_values, label="other", bottom=bottom, color="lightgray")

# Customize the plot
ax.set_ylabel('CPU usage %')
ax.set_xlabel('Thread number')
# ax.set_title('CPU Usage of Functions in Multi-threading', pad=20)

# Add margins to avoid label clipping
plt.margins(x=0.05, y=0.05)

# Define legend entries
legend_elements = []
for idx, function in enumerate(functions_to_be_printed):
    legend_elements.append(Patch(facecolor=colors[idx], label=function))

# Insert a separator after 'process_chrom'
separator = Patch(facecolor='white', edgecolor='white', label="───────────")
idx_process_chrom = functions_to_be_printed.index("process_chrom")
legend_elements.insert(idx_process_chrom + 1, separator)

# Add the "Other" category if needed
if np.any(other_values > 0):
    legend_elements.append(Patch(facecolor="lightgray", label="other"))

# Display legend
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.1, 1), borderaxespad=0., fontsize=12)

# Adjust layout to prevent overlap, increase padding around the plot
plt.tight_layout(pad=2.0)
plt.savefig('gprof.lcpan.threads.eps', format='eps')
