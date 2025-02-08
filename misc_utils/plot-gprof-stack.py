import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

data_t1 = """
 46.49      9.14     9.14       24     0.38     0.38  process_chrom
 15.87     12.26     3.12 15193490     0.00     0.00  find_boundaries
 12.51     14.72     2.46   580378     0.00     0.00  variate_sv
  4.48     15.60     0.88        1     0.88     1.25  print_ref_seq
  3.69     16.32     0.72 154754012     0.00     0.00  print_seq
  3.48     17.01     0.69 185136858     0.00     0.00  print_link
  2.87     17.57     0.56        1     0.56     8.00  t_read_vcf
  2.11     17.99     0.41 12438786     0.00     0.00  line_queue_push
  1.68     18.32     0.33 12438853     0.00     0.00  line_queue_pop
  1.68     18.65     0.33                             _init
  1.42     18.93     0.28        1     0.28     0.69  read_vcf
  1.27     19.18     0.25 14611057     0.00     0.00  variate_snp
  1.17     19.41     0.23 15193490     0.00     0.00  variate
  1.07     19.62     0.21        1     0.21     9.35  read_fasta
  0.15     19.65     0.03        1     0.03     0.03  free_ref_seq
  0.05     19.66     0.01                             suffix_to_prefix_overlap
  0.00     19.66     0.00        3     0.00     0.00  validate_file
  0.00     19.66     0.00        2     0.00     0.00  tpool_wait
  0.00     19.66     0.00        1     0.00     0.00  free_opt_arg
  0.00     19.66     0.00        1     0.00     0.00  line_queue_init
  0.00     19.66     0.00        1     0.00     0.00  parse_opts
  0.00     19.66     0.00        1     0.00     0.00  summarize
  0.00     19.66     0.00        1     0.00     0.00  tpool_add_work
  0.00     19.66     0.00        1     0.00     0.00  tpool_create
  0.00     19.66     0.00        1     0.00     0.00  tpool_destroy
  0.00     19.66     0.00        1     0.00     0.00  tpool_work_create
  0.00     19.66     0.00        1     0.00     0.00  tpool_work_destroy
  0.00     19.66     0.00        1     0.00     0.00  tpool_work_get
""".strip().split("\n")

data_t2 = """
 44.83      9.15     9.15       24     0.38     0.38  process_chrom
 18.50     12.93     3.77 11946106     0.00     0.00  find_boundaries
  9.58     14.88     1.96   438812     0.00     0.00  variate_sv
  5.29     15.96     1.08        1     1.08     1.65  print_ref_seq
  5.05     16.99     1.03 147042892     0.00     0.00  print_link
  3.65     17.73     0.74 124724236     0.00     0.00  print_seq
  2.87     18.32     0.58        2     0.29     4.21  t_read_vcf
  1.69     18.66     0.34  9155950     0.00     0.00  line_queue_push
  1.62     19.00     0.33 11362206     0.00     0.00  variate
  1.47     19.30     0.30                             _init
  1.45     19.59     0.29 10986369     0.00     0.00  variate_snp
  1.42     19.88     0.29        1     0.29     0.64  read_vcf
  1.35     20.16     0.28  8947314     0.00     0.00  line_queue_pop
  1.08     20.38     0.22        1     0.22     9.37  read_fasta
  0.12     20.40     0.03        1     0.03     0.03  free_ref_seq
  0.02     20.41     0.01        1     0.01     0.01  line_queue_init
  0.02     20.41     0.01                             suffix_to_prefix_overlap
  0.00     20.41     0.00        3     0.00     0.00  validate_file
  0.00     20.41     0.00        2     0.00     0.00  tpool_add_work
  0.00     20.41     0.00        2     0.00     0.00  tpool_wait
  0.00     20.41     0.00        2     0.00     0.00  tpool_work_create
  0.00     20.41     0.00        2     0.00     0.00  tpool_work_destroy
  0.00     20.41     0.00        2     0.00     0.00  tpool_work_get
  0.00     20.41     0.00        1     0.00     0.00  free_opt_arg
  0.00     20.41     0.00        1     0.00     0.00  parse_opts
  0.00     20.41     0.00        1     0.00     0.00  summarize
  0.00     20.41     0.00        1     0.00     0.00  tpool_create
  0.00     20.41     0.00        1     0.00     0.00  tpool_destroy
""".strip().split("\n")

data_t4 = """
 63.96      9.05     9.05       24     0.38     0.38  process_chrom
  9.29     10.37     1.31  8659750     0.00     0.00  find_boundaries
  7.77     11.46     1.10        1     1.10     1.47  print_ref_seq
  5.51     12.24     0.78   234179     0.00     0.00  variate_sv
  3.32     12.71     0.47 100233133     0.00     0.00  print_link
  2.30     13.04     0.33 86379488     0.00     0.00  print_seq
  1.48     13.25     0.21                             _init
  1.45     13.46     0.20        4     0.05     0.77  t_read_vcf
  1.41     13.65     0.20        1     0.20     9.25  read_fasta
  1.06     13.80     0.15  6325315     0.00     0.00  variate_snp
  0.78     13.91     0.11  6821139     0.00     0.00  variate
  0.64     14.01     0.09        1     0.09     0.14  read_vcf
  0.60     14.09     0.09  4967170     0.00     0.00  line_queue_pop
  0.35     14.14     0.05  5203563     0.00     0.00  line_queue_push
  0.04     14.14     0.01        1     0.01     0.01  free_ref_seq
  0.04     14.15     0.01                             suffix_to_prefix_overlap
  0.00     14.15     0.00        4     0.00     0.00  tpool_add_work
  0.00     14.15     0.00        4     0.00     0.00  tpool_work_create
  0.00     14.15     0.00        3     0.00     0.00  tpool_work_get
  0.00     14.15     0.00        3     0.00     0.00  validate_file
  0.00     14.15     0.00        2     0.00     0.00  tpool_wait
  0.00     14.15     0.00        2     0.00     0.00  tpool_work_destroy
  0.00     14.15     0.00        1     0.00     0.00  free_opt_arg
  0.00     14.15     0.00        1     0.00     0.00  line_queue_init
  0.00     14.15     0.00        1     0.00     0.00  parse_opts
  0.00     14.15     0.00        1     0.00     0.00  summarize
  0.00     14.15     0.00        1     0.00     0.00  tpool_create
  0.00     14.15     0.00        1     0.00     0.00  tpool_destroy
""".strip().split("\n")

data_t8 = """
 79.20      8.95     8.95       24     0.37     0.37  process_chrom
  8.58      9.92     0.97        1     0.97     1.44  print_ref_seq
  3.14     10.28     0.35 56044581     0.00     0.00  print_seq
  3.01     10.62     0.34    15423     0.00     0.00  variate_sv
  2.21     10.87     0.25 56880772     0.00     0.00  print_link
  1.42     11.03     0.16   700025     0.00     0.00  find_boundaries
  1.15     11.15     0.13        1     0.13     9.08  read_fasta
  0.80     11.24     0.09                             _init
  0.18     11.27     0.02        6     0.00     0.11  t_read_vcf
  0.09     11.28     0.01   400408     0.00     0.00  variate
  0.09     11.29     0.01   321219     0.00     0.00  line_queue_push
  0.09     11.29     0.01        1     0.01     0.02  read_vcf
  0.04     11.30     0.01        1     0.01     0.01  free_ref_seq
  0.00     11.30     0.00   385414     0.00     0.00  variate_snp
  0.00     11.30     0.00   323466     0.00     0.00  line_queue_pop
  0.00     11.30     0.00        8     0.00     0.00  tpool_work_get
  0.00     11.30     0.00        7     0.00     0.00  tpool_add_work
  0.00     11.30     0.00        7     0.00     0.00  tpool_work_create
  0.00     11.30     0.00        4     0.00     0.00  tpool_work_destroy
  0.00     11.30     0.00        3     0.00     0.00  validate_file
  0.00     11.30     0.00        1     0.00     0.00  free_opt_arg
  0.00     11.30     0.00        1     0.00     0.00  line_queue_init
  0.00     11.30     0.00        1     0.00     0.00  parse_opts
  0.00     11.30     0.00        1     0.00     0.00  summarize
  0.00     11.30     0.00        1     0.00     0.00  tpool_create
  0.00     11.30     0.00        1     0.00     0.00  tpool_destroy
  0.00     11.30     0.00        1     0.00     0.00  tpool_wait
""".strip().split("\n")

data_t16 = """
 65.48      9.01     9.01       24     0.38     0.38  process_chrom
  7.70     10.07     1.06   909682     0.00     0.00  find_boundaries
  6.87     11.02     0.94        1     0.94     1.74  print_ref_seq
  6.69     11.94     0.92    25203     0.00     0.00  variate_sv
  4.22     12.52     0.58 53679178     0.00     0.00  print_seq
  3.02     12.93     0.41 54940678     0.00     0.00  print_link
  1.53     13.14     0.21        1     0.21     9.22  read_fasta
  1.24     13.31     0.17       14     0.01     0.18  t_read_vcf
  0.87     13.43     0.12                             _init
  0.58     13.51     0.08   616622     0.00     0.00  variate
  0.51     13.58     0.07        1     0.07     0.12  read_vcf
  0.40     13.63     0.06   485313     0.00     0.00  line_queue_pop
  0.36     13.69     0.05   622757     0.00     0.00  variate_snp
  0.33     13.73     0.04   484709     0.00     0.00  line_queue_push
  0.22     13.76     0.03        1     0.03     0.03  free_ref_seq
  0.00     13.76     0.00       16     0.00     0.00  tpool_add_work
  0.00     13.76     0.00       16     0.00     0.00  tpool_work_create
  0.00     13.76     0.00       13     0.00     0.00  tpool_work_get
  0.00     13.76     0.00        4     0.00     0.00  tpool_work_destroy
  0.00     13.76     0.00        3     0.00     0.00  validate_file
  0.00     13.76     0.00        1     0.00     0.00  free_opt_arg
  0.00     13.76     0.00        1     0.00     0.00  line_queue_init
  0.00     13.76     0.00        1     0.00     0.00  parse_opts
  0.00     13.76     0.00        1     0.00     0.00  summarize
  0.00     13.76     0.00        1     0.00     0.00  tpool_create
  0.00     13.76     0.00        1     0.00     0.00  tpool_destroy
  0.00     13.76     0.00        1     0.00     0.00  tpool_wait
""".strip().split("\n")

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
