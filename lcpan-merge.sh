#!/bin/bash

log_file=$1

if [ ! -f "$log_file" ]; then
    echo "Error: Log file '$log_file' not found."
    exit 1
fi

program_mode=$(sed -n '1p' "$log_file")
input_file=$(sed -n '2p' "$log_file")
file_count=$(sed -n '3p' "$log_file")

if [ "$program_mode" == "vg" ]; then
    if ! [[ "$file_count" =~ ^[0-9]+$ ]]; then
        echo "Error: File count '$file_count' is not a valid number."
        exit 1
    fi

    rm -f "$input_file"

    for i in $(seq 0 $((file_count))); do
        segment_file="${input_file}.s.${i}"

        if [ -f "$segment_file" ]; then
            if ! cat "$segment_file" >> "$input_file"; then
                echo "Error appending $segment_file to $input_file"
                exit 1
            fi
            rm "$segment_file"
        else
            echo "Warning: File '$segment_file' does not exist. Skipping."
        fi
    done

    for i in $(seq 0 $((file_count))); do
        link_file="${input_file}.l.${i}"

        if [ -f "$link_file" ]; then
            if ! cat "$link_file" >> "$input_file"; then
                echo "Error appending $link_file to $input_file"
                exit 1
            fi
            rm "$link_file"
        else
            echo "Warning: File '$link_file' does not exist. Skipping."
        fi
    done

    if ! cat "${input_file}.p" >> "$input_file"; then
        echo "Error appending ${input_file}.p to $input_file"
        exit 1
    fi
    rm "${input_file}.p"
elif [ "$program_mode" == "vgx" ]; then
    if [ ! -f "$input_file" ]; then
        echo "Error: Input file '$input_file' is invalid or does not exist."
        exit 1
    fi

    if ! [[ "$file_count" =~ ^[0-9]+$ ]]; then
        echo "Error: File count '$file_count' is not a valid number."
        exit 1
    fi

    for i in $(seq 1 $((file_count))); do
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
fi
