#!/bin/bash

PROGRAM="./lcpan -f data/fasta_files/a.fa -v data/vcf_files/a.vcf -o out/testaaaaaaaaaaaa.rgfa -l 4"
INTERVAL=1

# Check for required dependencies
if ! command -v bc >/dev/null 2>&1; then
    echo "Error: 'bc' is not installed. Please install it to use this script."
    exit 1
fi

# Run the program in the background
$PROGRAM &
PID=$!
sleep 1

# Initialize stats
CPU_TOTAL=0
RAM_TOTAL=0
SAMPLES=0
MAX_CPU=0
MAX_RAM=0

# Trap to handle termination
trap "echo 'Monitoring stopped'; exit" SIGINT SIGTERM

while ps -p $PID > /dev/null; do
    # Capture current stats
    CPU=$(ps -p $PID -o %cpu --no-headers | awk '{print $1+0}')
    RAM=$(ps -p $PID -o rss --no-headers | awk '{print $1+0}')

    # Validate values
    if [[ -z $CPU || -z $RAM ]]; then
        sleep $INTERVAL
        continue
    fi

    # Update totals and max values
    CPU_TOTAL=$(echo "$CPU_TOTAL + $CPU" | bc)
    RAM_TOTAL=$(echo "$RAM_TOTAL + $RAM" | bc)
    MAX_CPU=$(echo "$MAX_CPU > $CPU ? $MAX_CPU : $CPU" | bc)
    MAX_RAM=$(echo "$MAX_RAM > $RAM ? $MAX_RAM : $RAM" | bc)
    SAMPLES=$((SAMPLES + 1))

    sleep $INTERVAL
done

# Check to avoid division by zero
if [ "$SAMPLES" -eq 0 ]; then
    echo "No data collected. Program might have terminated too quickly."
    exit 1
fi

# Calculate averages
AVG_CPU=$(echo "scale=2; $CPU_TOTAL / $SAMPLES" | bc)
AVG_RAM=$(echo "scale=2; $RAM_TOTAL / $SAMPLES" | bc)

echo "Maximum CPU Usage: $MAX_CPU%"
echo "Average CPU Usage: $AVG_CPU%"
echo "Maximum RAM Usage: $MAX_RAM KB"
echo "Average RAM Usage: $AVG_RAM KB"

