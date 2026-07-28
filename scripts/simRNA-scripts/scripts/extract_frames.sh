#!/bin/bash

# More robust version with better error handling and logging

BASE_DIR="tests_temps"
SIMRNA_EXEC="../soft/SimRNA_64bitIntel_Linux/SimRNA_trafl2pdbs"
LOG_FILE="simrna_processing.log"

# Initialize log file
echo "SimRNA processing started at $(date)" > "$LOG_FILE"

# Validation
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: Directory $BASE_DIR not found!" | tee -a "$LOG_FILE"
    exit 1
fi

if [ ! -f "$SIMRNA_EXEC" ]; then
    echo "Error: SimRNA executable not found at $SIMRNA_EXEC" | tee -a "$LOG_FILE"
    exit 1
fi

successful=0
failed=0
skipped=0

for subdir in "$BASE_DIR"/*/; do
    if [ -d "$subdir" ]; then
        dirname=$(basename "$subdir")
        echo "Processing: $dirname" | tee -a "$LOG_FILE"
        
        cd "$subdir"
        
        if [ -f "3WJ_0_1_3.trafl" ]; then
            if "$SIMRNA_EXEC" 3WJ_0_1_3-000001.pdb 3WJ_0_1_3.trafl 2>>"../$LOG_FILE"; then
                echo "  ✓ Success: $dirname" | tee -a "../$LOG_FILE"
                ((successful++))
            else
                echo "  ✗ Failed: $dirname" | tee -a "../$LOG_FILE"
                ((failed++))
            fi
        else
            echo "  ⚠ Skipped: $dirname (no .trafl file)" | tee -a "../$LOG_FILE"
            ((skipped++))
        fi
        
        cd - > /dev/null
    fi
done

echo "Summary: $successful successful, $failed failed, $skipped skipped" | tee -a "$LOG_FILE"
echo "Processing completed at $(date)" >> "$LOG_FILE"
