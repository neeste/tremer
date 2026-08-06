#!/usr/bin/env python3
import sys
import re
import os

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 tremer2sapp.py <input.txt> <output.txt>")
        sys.exit(1)
        
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    
    if not os.path.exists(in_file):
        print(f"Error: Could not find input file '{in_file}'")
        sys.exit(1)
        
    with open(in_file, 'r') as f:
        lines = f.readlines()
        
    snp_data = []
    str_data = []
    
    in_snp_data = False
    
    for line in lines:
        stripped = line.strip()
        
        if stripped.startswith('/SNPDATA'):
            in_snp_data = True
            continue
        elif stripped.startswith('/') and in_snp_data:
            in_snp_data = False
            
        if in_snp_data and stripped and not stripped.startswith('#'):
            snp_data.append(stripped)
            continue
            
        # Parse TREMER kit data lines
        # E.g., "133480 13 24 14 11 ..."
        parts = stripped.split()
        if len(parts) > 2:
            # Check if parts[1] is a number (a STR value)
            if parts[1].isdigit() or parts[1] == 'N' or '-' in parts[1]:
                # Looks like a kit data line
                str_data.append(stripped)
                
    with open(out_file, 'w') as f:
        f.write("/STRDATA\n")
        for st in str_data:
            f.write(f"{st}\n")
            
        if snp_data:
            f.write("\n/SNPDATA\n")
            for snp in snp_data:
                f.write(f"{snp}\n")
                
    print(f"Successfully converted TREMER file into SAPP format.")
    print(f"Exported {len(str_data)} kits.")

if __name__ == '__main__':
    main()
