import numpy as np

gfa_sizes_b = [, , , ]  # B
rgfa_sizes_b = [, , , ]  # B
ram_usage_kb = [, , , ]  # KB
execution_time_str = ["", "", "", ""]  # H:MM:SS.ss or MM:SS.ss
merge_time_str = ["", "", "", ""]  # H:MM:SS.ss or MM:SS.ss

# Convert RAM usage to GB
ram_usage_gb = [x / (1024 ** 2) for x in ram_usage_kb]  # 1 GB = 1024^2 KB

gfa_sizes_b = [x / (1024 ** 3) for x in gfa_sizes_b]
rgfa_sizes_b = [x / (1024 ** 3) for x in rgfa_sizes_b]

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

print(gfa_sizes_b)
print(rgfa_sizes_b)
print(ram_usage_gb)
print(execution_time_sec)