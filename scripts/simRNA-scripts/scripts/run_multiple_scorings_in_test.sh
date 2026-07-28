#!/bin/bash

# Enable job control explicitly
set -m

# Check if suffix argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <suffix>"
    echo "Example: $0 _ss005"
    echo "Example: $0 _ss010"
    exit 1
fi

# Get suffix from command line argument
suffix="$1"

# Create array to track PIDs
pids=()

for dir in *${suffix}; do
    if [ -d "$dir" ]; then
        # Extract base name
        base=$(basename "$dir" | grep -o 'C2_.*_S4_[0-9]\+')
        str_file="$dir/${base}.str"
        
        if [ -f "$str_file" ]; then
            echo "Starting job for $dir with base $base"
            
            # Start background job and capture PID
            (
                cd "$dir"
                python ../../../scripts/Analyze_scorings.py \
                    --name "$base" \
                    --min-temp 1 \
                    --max-temp 16 \
                    --structure "${base}.str" \
                    --summary ../medians_with_only_dots.csv \
                    --label "$base" \
                    --median \
                    --quiet
            ) &
            
            # Store the PID
            pids+=($!)
            echo "Started job with PID $!"
        else
            echo "Warning: $str_file not found, skipping."
        fi
    fi
done

echo "Started ${#pids[@]} jobs with PIDs: ${pids[*]}"

# Wait for all jobs to complete
for pid in "${pids[@]}"; do
    echo "Waiting for PID $pid"
    wait "$pid"
    echo "PID $pid finished"
done

echo "All jobs finished."