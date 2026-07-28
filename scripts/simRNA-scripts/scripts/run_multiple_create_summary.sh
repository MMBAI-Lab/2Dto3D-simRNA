#!/bin/bash
# filepath: /sdb-disk/pspina/RNA/simRNA/scripts/run_multiple_create_summary.sh

# Enable job control explicitly
set -m

# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <suffix> <summary_file> [--median|--mean] [--threshold VALUE]"
    echo "Example: $0 _ss005 results.csv --median"
    echo "Example: $0 _ss010 summary.csv --mean"
    echo "Example: $0 _ss005 results.csv --threshold 0.7"
    echo "Example: $0 _ss005 results.csv"
    exit 1
fi

# Get arguments
suffix="$1"
summary_file="$2"
shift 2  # Remove first two arguments, rest are optional flags

# Parse optional arguments
extra_args=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --median)
            extra_args="$extra_args --median"
            shift
            ;;
        --mean)
            extra_args="$extra_args --mean"
            shift
            ;;
        --threshold)
            extra_args="$extra_args --threshold $2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create array to track PIDs
pids=()

echo "Starting create_summary jobs for directories ending with '$suffix'"
echo "Summary file: $summary_file"
echo "Extra arguments: $extra_args"

for dir in *${suffix}; do
    if [ -d "$dir" ]; then
        # Extract base name
        base=$(basename "$dir" | grep -o 'C2_.*_S4_[0-9]\+')
        
        # Check if any accuracy CSV files exist in this directory
        if ls "$dir"/accuracy_TEMP*_pairs.csv 1> /dev/null 2>&1; then
            echo "Starting create_summary job for $dir with base $base"
            
            # Start background job and capture PID
            (
                cd "$dir"
                python ../../../scripts/create_summary.py \
                    --summary "../$summary_file" \
                    --label "$base" \
                    $extra_args
            ) &
            
            # Store the PID
            pids+=($!)
            echo "Started create_summary job with PID $!"
        else
            echo "Warning: No accuracy CSV files found in $dir, skipping."
        fi
    fi
done

echo "Started ${#pids[@]} create_summary jobs with PIDs: ${pids[*]}"

# Wait for all jobs to complete
for pid in "${pids[@]}"; do
    echo "Waiting for PID $pid"
    wait "$pid"
    echo "PID $pid finished"
done

echo "All create_summary jobs finished."
echo "Summary file '$summary_file' should now contain all results."