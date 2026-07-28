#!/usr/bin/env python3
"""
SimRNA Temperature Sweep Setup Script

Creates a folder structure for SimRNA temperature sweep simulations.
Each temperature folder contains:
- Sequence file (like 3WJ_0_1_3)
- Structure file (like 3WJ_0_1_3.str) 
- Configuration file (sim_config.in) with temperature-specific settings
- Symlink to SimRNA data directory

Example usage:
  python Temp_sweep.py --output-dir temp_sweep_test \
                       --sequence-file sequences/3WJ_0_1_3 \
                       --structure-file structures/3WJ_0_1_3.str \
                       --config-template templates/sim_config_template.in \
                       --temp-start 1.20 --temp-end 2.00 --temp-step 0.05 \
                       --create-run-script

  # For a quick test with fewer temperatures:
  python Temp_sweep.py --output-dir quick_test \
                       --sequence-file my_sequence.txt \
                       --structure-file my_sequence.str \
                       --config-template config_template.in \
                       --temp-start 1.40 --temp-end 1.60 --temp-step 0.10
"""

import os
import argparse
import shutil
import subprocess
from pathlib import Path

def create_temp_folders(output_dir, sequence_file, structure_file, config_template, 
                       temp_start=1.20, temp_end=2.00, temp_step=0.05, 
                       simrna_data_path="../../soft/SimRNA_64bitIntel_Linux/data"):
    """
    Create temperature sweep folder structure for SimRNA.
    
    Args:
        output_dir: Base directory for temperature folders
        sequence_file: Path to sequence file (like 3WJ_0_1_3)
        structure_file: Path to structure file (like 3WJ_0_1_3.str)
        config_template: Path to template sim_config.in file
        temp_start: Starting temperature (default: 1.20)
        temp_end: Ending temperature (default: 2.00)
        temp_step: Temperature step (default: 0.05)
        simrna_data_path: Relative path to SimRNA data directory
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get base names for files
    seq_basename = os.path.basename(sequence_file)
    struct_basename = os.path.basename(structure_file)
    
    # Read template config file
    with open(config_template, 'r') as f:
        config_template_content = f.read()
    
    # Generate temperature range
    current_temp = temp_start
    temp_count = 0
    
    print(f"Creating temperature sweep from {temp_start} to {temp_end} with step {temp_step}")
    print(f"Output directory: {output_dir}")
    print(f"SimRNA data symlink: {simrna_data_path}")
    print()
    
    while current_temp <= temp_end + 1e-10:  # Small epsilon for floating point comparison
        temp_folder = f"temp_{current_temp:.2f}"
        temp_dir = os.path.join(output_dir, temp_folder)
        
        print(f"Creating {temp_folder}...")
        
        # Create temperature directory
        os.makedirs(temp_dir, exist_ok=True)
        
        # Copy sequence file
        dest_seq = os.path.join(temp_dir, seq_basename)
        shutil.copy2(sequence_file, dest_seq)
        
        # Copy structure file
        dest_struct = os.path.join(temp_dir, struct_basename)
        shutil.copy2(structure_file, dest_struct)
        
        # Create temperature-specific config file
        config_content = config_template_content.replace("INIT_TEMP", f"INIT_TEMP {current_temp:.2f}")
        config_content = config_content.replace("FINAL_TEMP", f"FINAL_TEMP {current_temp:.2f}")
        
        # Handle case where template might already have specific temperature values
        lines = config_content.split('\n')
        new_lines = []
        init_temp_set = False
        final_temp_set = False
        
        for line in lines:
            if line.strip().startswith('INIT_TEMP'):
                new_lines.append(f"INIT_TEMP {current_temp:.2f}")
                init_temp_set = True
            elif line.strip().startswith('FINAL_TEMP'):
                new_lines.append(f"FINAL_TEMP {current_temp:.2f}")
                final_temp_set = True
            else:
                new_lines.append(line)
        
        # If INIT_TEMP or FINAL_TEMP weren't found, add them
        if not init_temp_set:
            new_lines.insert(-1, f"INIT_TEMP {current_temp:.2f}")
        if not final_temp_set:
            new_lines.insert(-1, f"FINAL_TEMP {current_temp:.2f}")
        
        config_content = '\n'.join(new_lines)
        
        # Write config file
        config_path = os.path.join(temp_dir, "sim_config.in")
        with open(config_path, 'w') as f:
            f.write(config_content)
        
        # Create symlink to data directory
        data_link = os.path.join(temp_dir, "data")
        if os.path.exists(data_link) or os.path.islink(data_link):
            os.remove(data_link)  # Remove existing symlink
        
        try:
            # Create relative symlink
            os.symlink(simrna_data_path, data_link)
            print(f"  ✓ Created symlink: data -> {simrna_data_path}")
        except OSError as e:
            print(f"  ✗ Failed to create symlink: {e}")
        
        temp_count += 1
        current_temp += temp_step
    
    print(f"\nCompleted! Created {temp_count} temperature folders.")
    print(f"Temperature range: {temp_start:.2f} to {current_temp - temp_step:.2f}")
    
    # Create a summary file
    summary_file = os.path.join(output_dir, "temperature_sweep_info.txt")
    with open(summary_file, 'w') as f:
        f.write(f"SimRNA Temperature Sweep Configuration\n")
        f.write(f"======================================\n\n")
        f.write(f"Output directory: {os.path.abspath(output_dir)}\n")
        f.write(f"Temperature range: {temp_start:.2f} to {current_temp - temp_step:.2f}\n")
        f.write(f"Temperature step: {temp_step:.2f}\n")
        f.write(f"Number of folders: {temp_count}\n")
        f.write(f"Sequence file: {seq_basename}\n")
        f.write(f"Structure file: {struct_basename}\n")
        f.write(f"Config template: {os.path.basename(config_template)}\n")
        f.write(f"SimRNA data path: {simrna_data_path}\n")
        f.write(f"\nFolders created:\n")
        
        current_temp = temp_start
        while current_temp <= temp_end + 1e-10:
            f.write(f"  temp_{current_temp:.2f}/\n")
            current_temp += temp_step
    
    print(f"Summary written to: {summary_file}")

def create_run_script(output_dir, simrna_executable="../soft/SimRNA_64bitIntel_Linux/SimRNA"):
    """Create a script to run SimRNA in all temperature directories."""
    
    script_path = os.path.join(output_dir, "run_all_temps.sh")
    
    with open(script_path, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Auto-generated script to run SimRNA for all temperatures\n\n")
        f.write("echo \"Starting SimRNA temperature sweep at $(date)\"\n\n")
        
        # Find all temperature directories
        f.write("for temp_dir in temp_*/; do\n")
        f.write("    if [ -d \"$temp_dir\" ]; then\n")
        f.write("        echo \"Running SimRNA in $temp_dir\"\n")
        f.write("        cd \"$temp_dir\"\n")
        f.write("        \n")
        f.write("        # Get the base name of sequence file (without path)\n")
        f.write("        SEQ_FILE=$(ls | grep -v '.str$' | grep -v 'sim_config.in' | grep -v 'data' | head -1)\n")
        f.write("        STRUCT_FILE=\"${SEQ_FILE}.str\"\n")
        f.write("        \n")
        f.write("        echo \"  Sequence: $SEQ_FILE\"\n")
        f.write("        echo \"  Structure: $STRUCT_FILE\"\n")
        f.write("        echo \"  Output: $SEQ_FILE\"\n")
        f.write("        \n")
        f.write(f"        {simrna_executable} -s $SEQ_FILE -S $STRUCT_FILE -c sim_config.in -o $SEQ_FILE >& ${{SEQ_FILE}}.log &\n")
        f.write("        \n")
        f.write("        cd ..\n")
        f.write("        echo \"Started SimRNA for $temp_dir (running in background)\"\n")
        f.write("        echo \"  Log file: $temp_dir/${SEQ_FILE}.log\"\n")
        f.write("    fi\n")
        f.write("done\n\n")
        f.write("echo \"All SimRNA jobs started at $(date)\"\n")
        f.write("echo \"Monitor progress with: tail -f temp_*/[sequence_name].log\"\n")
        f.write("echo \"Check running jobs with: ps aux | grep SimRNA\"\n")
    
    os.chmod(script_path, 0o755)
    print(f"Created run script: {script_path}")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Create SimRNA temperature sweep folder structure",
        epilog="""
Example commands:

  # Full temperature sweep (1.20 to 2.00, step 0.05):
  python %(prog)s --output-dir full_temp_sweep \\
                  --sequence-file sequences/3WJ_0_1_3 \\
                  --structure-file sequences/3WJ_0_1_3.str \\
                  --config-template templates/sim_config.in \\
                  --create-run-script

  # Quick test with fewer temperatures:
  python %(prog)s --output-dir quick_test \\
                  --sequence-file my_rna.seq \\
                  --structure-file my_rna.str \\
                  --config-template config.in \\
                  --temp-start 1.40 --temp-end 1.60 --temp-step 0.10

  # Custom temperature range:
  python %(prog)s --output-dir custom_sweep \\
                  --sequence-file hairpin.seq \\
                  --structure-file hairpin.str \\
                  --config-template sim_config_template.in \\
                  --temp-start 1.50 --temp-end 1.80 --temp-step 0.02 \\
                  --simrna-data-path /path/to/SimRNA/data \\
                  --create-run-script

After running, use the generated run script:
  cd output_directory && ./run_all_temps.sh
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--output-dir", required=True,
                       help="Output directory for temperature folders")
    parser.add_argument("--sequence-file", required=True,
                       help="Path to sequence file (like 3WJ_0_1_3)")
    parser.add_argument("--structure-file", required=True,
                       help="Path to structure file (like 3WJ_0_1_3.str)")
    parser.add_argument("--config-template", required=True,
                       help="Path to template sim_config.in file")
    
    parser.add_argument("--temp-start", type=float, default=1.20,
                       help="Starting temperature (default: 1.20)")
    parser.add_argument("--temp-end", type=float, default=2.00,
                       help="Ending temperature (default: 2.00)")
    parser.add_argument("--temp-step", type=float, default=0.05,
                       help="Temperature step (default: 0.05)")
    
    parser.add_argument("--simrna-data-path", 
                       default="../../soft/SimRNA_64bitIntel_Linux/data",
                       help="Relative path to SimRNA data directory")
    parser.add_argument("--simrna-executable",
                       default="../../soft/SimRNA_64bitIntel_Linux/SimRNA",
                       help="Path to SimRNA executable for run script")
    
    parser.add_argument("--create-run-script", action="store_true",
                       help="Also create a script to run all simulations")
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Validate input files
    for file_path, name in [(args.sequence_file, "sequence file"), 
                           (args.structure_file, "structure file"),
                           (args.config_template, "config template")]:
        if not os.path.exists(file_path):
            print(f"Error: {name} not found: {file_path}")
            return 1
    
    # Create temperature folders
    create_temp_folders(
        output_dir=args.output_dir,
        sequence_file=args.sequence_file,
        structure_file=args.structure_file,
        config_template=args.config_template,
        temp_start=args.temp_start,
        temp_end=args.temp_end,
        temp_step=args.temp_step,
        simrna_data_path=args.simrna_data_path
    )
    
    # Optionally create run script
    if args.create_run_script:
        create_run_script(args.output_dir, args.simrna_executable)
    
    print("\n" + "="*60)
    print("NEXT STEPS:")
    print("="*60)
    print("\n1. To run individual simulations:")
    print(f"   cd {args.output_dir}/temp_1.20")
    print("   ../../soft/SimRNA_64bitIntel_Linux/SimRNA -s [sequence_file] -S [structure_file] -c sim_config.in")
    
    if args.create_run_script:
        print("\n2. To run all simulations at once:")
        print(f"   cd {args.output_dir}")
        print("   ./run_all_temps.sh")
        print("\n3. Monitor progress:")
        print("   tail -f temp_*/[sequence_name].log")
        print("   ps aux | grep SimRNA")
    
    print(f"\n4. Temperature folders created: {args.temp_start:.2f} to {args.temp_end:.2f} (step {args.temp_step:.2f})")
    print(f"   Total folders: {int((args.temp_end - args.temp_start) / args.temp_step) + 1}")

if __name__ == "__main__":
    main()