#!/usr/bin/env python3
import sys
import os
import glob

def analyze_file_structure(filename):
    """
    Analyze the file structure in detail to identify specific issues
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Keep original lines for analysis
        original_lines = [line.rstrip('\n\r') for line in lines]
        # Also get stripped lines for processing
        stripped_lines = [line.strip() for line in lines if line.strip()]
        
        print(f"\n📊 Detailed analysis of {filename}:")
        print(f"   Total lines (including empty): {len(original_lines)}")
        print(f"   Non-empty lines: {len(stripped_lines)}")
        print(f"   Empty lines: {len(original_lines) - len(stripped_lines)}")
        
        # Show first few lines for context
        print(f"\n📝 First 10 lines preview:")
        for i, line in enumerate(original_lines[:10]):
            line_type = "EMPTY" if not line.strip() else f"{len(line.split())} numbers"
            print(f"   Line {i+1:3d}: [{line_type:10s}] '{line[:50]}{'...' if len(line) > 50 else ''}'")
        
        if len(stripped_lines) == 0:
            return False, "File contains no non-empty lines"
        
        # Analyze pattern
        print(f"\n🔍 Pattern analysis:")
        line_lengths = [len(line.split()) for line in stripped_lines]
        unique_lengths = list(set(line_lengths))
        unique_lengths.sort()
        
        print(f"   Unique line lengths (number of values): {unique_lengths}")
        
        # Count occurrences of each length
        length_counts = {}
        for length in line_lengths:
            length_counts[length] = length_counts.get(length, 0) + 1
        
        print(f"   Length distribution:")
        for length in sorted(length_counts.keys()):
            count = length_counts[length]
            percentage = (count / len(stripped_lines)) * 100
            print(f"     {length:3d} numbers: {count:5d} lines ({percentage:5.1f}%)")
        
        # Check if it follows expected pattern
        return analyze_pattern(stripped_lines, line_lengths)
        
    except Exception as e:
        return False, f"Error reading file: {e}"

def analyze_pattern(lines, line_lengths):
    """
    Analyze if the file follows the expected alternating pattern
    """
    print(f"\n🎯 Pattern validation:")
    
    # Expected pattern: alternating lines of 5 numbers and many numbers
    if len(lines) % 2 != 0:
        odd_line_num = len(lines)
        print(f"   ❌ Odd number of lines ({len(lines)})")
        print(f"   💡 Last line {odd_line_num} might be incomplete or extra")
        print(f"   💡 Content of last line: '{lines[-1][:100]}{'...' if len(lines[-1]) > 100 else ''}'")
        
        # Check if removing last line would fix the pattern
        if analyze_alternating_pattern(lines[:-1], line_lengths[:-1]):
            return False, f"Odd number of lines ({len(lines)}). Removing last line would fix the format."
        else:
            return False, f"Odd number of lines ({len(lines)}). Pattern issues beyond just the extra line."
    
    return analyze_alternating_pattern(lines, line_lengths)

def analyze_alternating_pattern(lines, line_lengths):
    """
    Check if lines follow the expected alternating pattern
    """
    num_entries = len(lines) // 2
    expected_second_line_length = None
    issues = []
    
    print(f"   Checking {num_entries} entries for alternating pattern...")
    
    for i in range(0, len(lines), 2):
        entry_num = (i // 2) + 1
        first_line_length = line_lengths[i]
        second_line_length = line_lengths[i + 1]
        
        # Check first line (should have 5 numbers)
        if first_line_length != 5:
            issues.append(f"Entry {entry_num}: First line has {first_line_length} numbers (expected 5)")
            if len(issues) <= 5:  # Show first few issues
                print(f"   ❌ Entry {entry_num}: First line has {first_line_length} numbers, expected 5")
                print(f"      Content: '{lines[i][:80]}{'...' if len(lines[i]) > 80 else ''}'")
        
        # Check second line consistency
        if expected_second_line_length is None:
            expected_second_line_length = second_line_length
            print(f"   ℹ️  Expected second line length set to: {expected_second_line_length}")
        elif second_line_length != expected_second_line_length:
            issues.append(f"Entry {entry_num}: Second line has {second_line_length} numbers (expected {expected_second_line_length})")
            if len(issues) <= 5:  # Show first few issues
                print(f"   ❌ Entry {entry_num}: Second line has {second_line_length} numbers, expected {expected_second_line_length}")
        
        # Validate that lines contain only numbers
        try:
            [float(x) for x in lines[i].split()]
            [float(x) for x in lines[i + 1].split()]
        except ValueError as e:
            issues.append(f"Entry {entry_num}: Non-numeric values found")
            if len(issues) <= 5:
                print(f"   ❌ Entry {entry_num}: Non-numeric values: {e}")
    
    if len(issues) > 5:
        print(f"   ... and {len(issues) - 5} more issues")
    
    if issues:
        return False, f"Found {len(issues)} pattern issues"
    else:
        print(f"   ✅ All entries follow expected pattern")
        return True, f"Format OK: {num_entries} entries, first lines have 5 numbers, second lines have {expected_second_line_length} numbers each"

def check_file_format(filename):
    """
    Enhanced format checker with detailed analysis
    """
    is_valid, message = analyze_file_structure(filename)
    return is_valid, message

def main():
    if len(sys.argv) > 1:
        # Check specific files provided as arguments
        files = sys.argv[1:]
    else:
        # Auto-detect level*_temp files in current directory
        files = glob.glob("level*_temp")
        if not files:
            print("No level*_temp files found in current directory")
            print("Usage: python check_format.py [file1] [file2] ...")
            sys.exit(1)
    
    print(f"Checking {len(files)} file(s):")
    print("=" * 80)
    
    all_valid = True
    
    for filename in sorted(files):
        if not os.path.exists(filename):
            print(f"❌ {filename}: File not found")
            all_valid = False
            continue
        
        is_valid, message = check_file_format(filename)
        
        if is_valid:
            print(f"\n✅ {filename}: {message}")
        else:
            print(f"\n❌ {filename}: {message}")
            all_valid = False
        
        print("-" * 80)
    
    if all_valid:
        print("\n🎉 All files have correct format!")
    else:
        print("\n⚠️  Some files have format issues")
        sys.exit(1)

if __name__ == "__main__":
    main()