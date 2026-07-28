#!/bin/bash
# filepath: /sdb-disk/pspina/RNA/simRNA/process_trajectories.sh

# Combined script to extract frames and create trajectory files for all subdirectories
# Usage: ./process_trajectories.sh <base_filename>
# Example: ./process_trajectories.sh 3WJ_BL

# Check if base filename is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <base_filename>"
    echo "Example: $0 3WJ_BL"
    exit 1
fi

BASE_NAME="$1"
SCRIPT_DIR="$(dirname "$0")"
SIMRNA_EXEC="$SCRIPT_DIR/../soft/SimRNA_64bitIntel_Linux/SimRNA_trafl2pdbs"
LOG_FILE="trajectory_processing.log"

# Initialize log file
echo "Trajectory processing started at $(date)" > "$LOG_FILE"
echo "Processing directory: $(pwd)" | tee -a "$LOG_FILE"
echo "Base filename: $BASE_NAME" | tee -a "$LOG_FILE"

create_trajectory_pdb() {
    local output_file="trajectory.pdb"
    local model_num=1
    
    # Find all PDB files with the pattern, sort them numerically
    local pdb_files=($(ls ${BASE_NAME}-*.pdb 2>/dev/null | sort -V))
    
    # Check if any PDB files exist
    if [ ${#pdb_files[@]} -eq 0 ]; then
        echo "No PDB files found matching pattern ${BASE_NAME}-*.pdb"
        return 1
    fi
    
    echo "Found ${#pdb_files[@]} PDB frames"
    
    # Create the trajectory file
    > "$output_file"  # Clear the file
    
    for pdb_file in "${pdb_files[@]}"; do
        if [ -f "$pdb_file" ]; then
            echo "MODEL     $(printf "%4d" $model_num)" >> "$output_file"
            
            # Copy ATOM/HETATM records, preserving exact formatting
            awk '/^ATOM|^HETATM/ {print}' "$pdb_file" >> "$output_file"
            
            echo "ENDMDL" >> "$output_file"
            ((model_num++))
        fi
    done
    
    echo "END" >> "$output_file"
    echo "Created trajectory.pdb with $((model_num-1)) frames"
    
    # Verify the trajectory was created successfully
    if [ -s "$output_file" ]; then
        local frame_count=$(grep -c "^MODEL" "$output_file")
        echo "Trajectory verification: $frame_count models in trajectory.pdb"
        return 0
    else
        echo "Error: trajectory.pdb is empty or was not created"
        return 1
    fi
}

successful=0
failed=0
skipped=0

echo "Processing subdirectories in current directory..." | tee -a "$LOG_FILE"

for subdir in */; do
    if [ -d "$subdir" ]; then
        dirname=$(basename "$subdir")
        echo "Processing: $dirname" | tee -a "$LOG_FILE"
        
        cd "$subdir"
        
        # Step 1: Extract frames if .trafl file exists
        if [ -f "${BASE_NAME}.trafl" ]; then
            echo "  Extracting frames..." | tee -a "../$LOG_FILE"
            
            # Check if reference structure exists
            if [ ! -f "${BASE_NAME}-000001.pdb" ]; then
                echo "  ✗ Reference structure ${BASE_NAME}-000001.pdb not found" | tee -a "../$LOG_FILE"
                ((failed++))
                cd - > /dev/null
                continue
            fi
            
            # Extract frames from trajectory
            if "$SIMRNA_EXEC" ${BASE_NAME}-000001.pdb ${BASE_NAME}.trafl : 2>>"../$LOG_FILE"; then
                echo "  ✓ Frame extraction successful" | tee -a "../$LOG_FILE"
                
                # Count extracted frames
                frame_count=$(ls ${BASE_NAME}-*.pdb 2>/dev/null | wc -l)
                echo "  Extracted $frame_count frames" | tee -a "../$LOG_FILE"
                
                # Step 2: Create multi-model PDB
                echo "  Creating trajectory PDB..." | tee -a "../$LOG_FILE"
                if create_trajectory_pdb >> "../$LOG_FILE" 2>&1; then
                    echo "  ✓ Trajectory PDB created" | tee -a "../$LOG_FILE"
                    
                    # Step 3: Clean up frame files but keep the reference structure
                    echo "  Cleaning up extracted frame files..." | tee -a "../$LOG_FILE"
                    
                    # Remove all frames except the first one (reference structure)
                    for frame_file in ${BASE_NAME}-*.pdb; do
                        if [ "$frame_file" != "${BASE_NAME}-000001.pdb" ]; then
                            rm -f "$frame_file"
                        fi
                    done
                    
                    # Clean up .ss_detected files
                    echo "  Cleaning up secondary structure files..." | tee -a "../$LOG_FILE"
                    #rm -f ${BASE_NAME}-*.ss_detected
                    
                    echo "  ✓ Success: $dirname (trajectory.pdb ready, reference structure preserved)" | tee -a "../$LOG_FILE"
                    ((successful++))
                else
                    echo "  ✗ Failed to create trajectory PDB: $dirname" | tee -a "../$LOG_FILE"
                    ((failed++))
                fi
            else
                echo "  ✗ Frame extraction failed: $dirname" | tee -a "../$LOG_FILE"
                ((failed++))
            fi
        else
            echo "  ⚠ Skipped: $dirname (no ${BASE_NAME}.trafl file)" | tee -a "../$LOG_FILE"
            ((skipped++))
        fi
        
        cd - > /dev/null
    fi
done

echo "" | tee -a "$LOG_FILE"
echo "Summary:" | tee -a "$LOG_FILE"
echo "  $successful successful (trajectory.pdb created)" | tee -a "$LOG_FILE"
echo "  $failed failed" | tee -a "$LOG_FILE"
echo "  $skipped skipped (no ${BASE_NAME}.trafl file)" | tee -a "$LOG_FILE"
echo "Processing completed at $(date)" | tee -a "$LOG_FILE"

# Summary of what was created
if [ $successful -gt 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "Created trajectory files:" | tee -a "$LOG_FILE"
    find . -name "trajectory.pdb" -exec ls -lh {} \; | tee -a "$LOG_FILE"
fi

# Provide usage information for the trajectory files
echo "" | tee -a "$LOG_FILE"
echo "Usage information:" | tee -a "$LOG_FILE"
echo "  - trajectory.pdb files contain all frames as MODEL/ENDMDL blocks" | tee -a "$LOG_FILE"
echo "  - Reference structure (${BASE_NAME}-000001.pdb) preserved in each directory" | tee -a "$LOG_FILE"
echo "  - Load trajectory.pdb in VMD, PyMOL, or other molecular viewers" | tee -a "$LOG_FILE"