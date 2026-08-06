#!/usr/bin/env python3
import sys
import re
import os

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 sapp2tremer.py <input.txt> <output.txt>")
        sys.exit(1)
        
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    
    if not os.path.exists(in_file):
        print(f"Error: Could not find input file '{in_file}'")
        sys.exit(1)
        
    kits = []
    
    with open(in_file, 'r') as f:
        content = f.read()
        
    # Extract kit IDs from /STRDATA block
    strdata_match = re.search(r'/STRDATA\s+(.*?)(?=\n/|$)', content, re.DOTALL | re.IGNORECASE)
    if strdata_match:
        lines = strdata_match.group(1).strip().split('\n')
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if parts:
                kits.append(parts[0])
                
    if not kits:
        print("Error: No kits found in /STRDATA block.")
        sys.exit(1)
        
    with open(out_file, 'w') as f:
        f.write("/GROUPS\n")
        f.write(f"SAPP_Import ({' '.join(kits)})\n\n")
        f.write(content)
        
    print(f"Successfully converted SAPP file into TREMER format.")
    print(f"Imported {len(kits)} kits into group 'SAPP_Import'.")

if __name__ == '__main__':
    main()
