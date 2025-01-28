#!/bin/bash

log_file="lcp_vg.log"

if [ ! -f "$log_file" ]; then
    echo "Error: Log file '$log_file' not found."
    exit 1
fi

input_file=$(sed -n '1p' "$log_file")
file_count=$(sed -n '2p' "$log_file")

if [ -z "$input_file" ] || [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' is invalid or does not exist."
    exit 1
fi

if ! [[ "$file_count" =~ ^[0-9]+$ ]]; then
    echo "Error: File count '$file_count' is not a valid number."
    exit 1
fi

for i in $(seq 0 $((file_count - 1))); do
    file="${input_file}.${i}"

    if [ -f "$file" ]; then
        if ! cat "$file" >> "$input_file"; then
            echo "Error appending $file to $input_file"
            exit 1
        fi
        rm "$file"
    else
        echo "Warning: File '$file' does not exist. Skipping."
    fi
done