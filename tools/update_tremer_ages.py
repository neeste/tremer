import urllib.request
import json
import sys
import time
import os
import glob
import ssl # Added to handle Mac PyInstaller SSL bugs

def detect_haplogroup(lines):
    """Scans the file to auto-detect the Haplogroup prefix (e.g., R, I, E)."""
    in_snptree = False
    for line in lines:
        raw_line = line.rstrip('\r\n')
        if '*' in raw_line:
            raw_line = raw_line[:raw_line.find('*')]
        stripped = raw_line.strip()
        
        if not stripped: continue
        if stripped.startswith('/'):
            in_snptree = (stripped == '/SNPTREE')
            continue
            
        if in_snptree:
            parts = stripped.split()
            if len(parts) >= 3:
                try:
                    int(parts[1]) 
                    child = parts[0]
                    parent = parts[2]
                    
                    for snp in [child, parent]:
                        if '-' in snp:
                            prefix = snp.split('-')[0]
                            if len(prefix) == 1 and prefix.isalpha():
                                return prefix.upper()
                                
                    for snp in [child, parent]:
                        if len(snp) == 1 and snp.isalpha():
                            return snp.upper()
                except ValueError:
                    pass
    return None 

def fetch_ftdna_age(snp_name, primary_haplo, retries=4):
    """Fetches the latest TMRCA mean age, specifically handling 429 Rate Limits."""
    clean_snp = snp_name
    if '-' in clean_snp:
        prefix = clean_snp.split('-')[0]
        if len(prefix) == 1 and prefix.isalpha():
            clean_snp = clean_snp.split('-')[1]
            
    prefixes_to_try = [f"{primary_haplo}-"] if primary_haplo else ['R-', 'I-', 'E-', 'J-', 'G-']
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
    }
    
    # Create an unverified SSL context to bypass PyInstaller certificate bugs
    ctx = ssl._create_unverified_context()
    
    for prefix in prefixes_to_try:
        url = f"https://discover.familytreedna.com/resources/y-dna/{prefix}{clean_snp}.json"
        last_error_code = None
        
        for attempt in range(retries):
            try:
                req = urllib.request.Request(url, headers=headers)
                # Pass the unverified context here!
                with urllib.request.urlopen(req, context=ctx) as response:
                    data = json.loads(response.read().decode('utf-8'))
                    if 'time' in data and data['time'] and 'tmrca' in data['time']:
                        return round(data['time']['tmrca'].get('mean', 0))
                    return 0
                    
            except urllib.error.HTTPError as e:
                last_error_code = e.code
                if e.code == 404:
                    break 
                else:
                    wait_time = 15 * (attempt + 1)
                    print(f"     [Server blocked (HTTP {e.code}). Pausing {wait_time}s to cool down...]")
                    time.sleep(wait_time)
            except Exception as e:
                wait_time = 15 * (attempt + 1)
                # We now print exactly what the exception is so we aren't guessing!
                print(f"     [Connection error: {e}. Pausing {wait_time}s to cool down...]")
                time.sleep(wait_time)
        
        if last_error_code != 404:
            return None 
            
    return -1

def process_tremer_file(filepath):
    print(f"--- Processing {filepath} ---")
    
    if not os.path.exists(filepath):
        print(f"File '{filepath}' not found!")
        return

    if "updated" in filepath.lower():
        print("Skipping file (appears to be an already-updated output).")
        return

    # PHASE 1: Quick read to collect all SNPs
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        lines = f.readlines()

    detected_haplo = detect_haplogroup(lines)
    
    snps_to_fetch = []
    in_snptree = False

    for line in lines:
        raw_line = line.rstrip('\r\n')
        if '*' in raw_line:
            data_part = raw_line[:raw_line.find('*')]
        else:
            data_part = raw_line
            
        stripped = data_part.strip()
        
        if stripped.startswith('/'):
            in_snptree = (stripped == '/SNPTREE')
            continue
            
        if in_snptree and stripped:
            parts = stripped.split()
            if len(parts) >= 3:
                child_snp = parts[0]
                middle_val = parts[1]
                
                try:
                    old_age = int(middle_val)
                    snps_to_fetch.append((child_snp, old_age))
                except ValueError:
                    pass

    if not snps_to_fetch:
        print("No SNPs found in /SNPTREE section.")
        return

    # PHASE 2: Fetch ages (SLOW, takes minutes)
    print(f"Found {len(snps_to_fetch)} SNPs. Fetching ages from FTDNA...")
    new_ages_dict = {}
    updates_made = 0
    
    for child_snp, old_age in snps_to_fetch:
        new_age = fetch_ftdna_age(child_snp, detected_haplo)
        
        if new_age == -1:
            print(f" [?] {child_snp} not found on FTDNA Discover. Keeping old age.")
        elif new_age is not None and new_age > 0:
            if new_age != old_age:
                print(f" [*] UPDATED {child_snp}: {old_age} -> {new_age}")
                new_ages_dict[child_snp] = new_age
                updates_made += 1
        else:
            print(f" [!] Error fetching {child_snp} (Rate Limit/Connection). Keeping old age.")
            
        time.sleep(2.0)

    if updates_made == 0:
        print("No new ages found. Skipping file write.")
        return

    # PHASE 3: Re-read the file (in case user edited it while fetching) and apply updates instantly
    print("Re-reading file to apply updates without overwriting recent manual edits...")
    
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        latest_lines = f.readlines()
        
    updated_lines = []
    in_snptree = False
    
    for line in latest_lines:
        raw_line = line.rstrip('\r\n')
        
        if '*' in raw_line:
            split_idx = raw_line.find('*')
            data_part = raw_line[:split_idx]
            comment_part = raw_line[split_idx:]
        else:
            data_part = raw_line
            comment_part = ""
            
        stripped = data_part.strip()
        
        if not stripped:
            updated_lines.append(line)
            continue
            
        if stripped.startswith('/'):
            in_snptree = (stripped == '/SNPTREE')
            updated_lines.append(line)
            continue
            
        if in_snptree and stripped:
            parts = stripped.split()
            if len(parts) >= 3:
                child_snp = parts[0]
                middle_val = parts[1]
                parent_snp = parts[2]
                
                try:
                    int(middle_val) # Verify it's a valid tree line
                    if child_snp in new_ages_dict:
                        new_age = new_ages_dict[child_snp]
                        new_data_str = f"{child_snp:<10} {new_age:<5} {parent_snp:<10}"
                        if comment_part:
                            updated_lines.append(f"{new_data_str} {comment_part}\n")
                        else:
                            updated_lines.append(f"{new_data_str}\n")
                        continue
                except ValueError:
                    pass
                    
        updated_lines.append(line)
        
    # Overwrite the actual input file safely
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(updated_lines)

    print(f"Done! Made {updates_made} age updates in {filepath}")

if __name__ == "__main__":
    files_to_process = []

    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            expanded_files = glob.glob(arg)
            if not expanded_files:
                print(f"Warning: No files matched '{arg}'")
            files_to_process.extend(expanded_files)
            
    else:
        print("No files provided. Opening File Picker...")
        
        import tkinter as tk
        from tkinter import filedialog
        
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        
        selected_files = filedialog.askopenfilenames(
            title="Select TREMER Input Files to Update",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
        )
        
        files_to_process = list(selected_files)

    files_to_process = list(set(files_to_process)) 
    
    if not files_to_process:
        print("\nNo valid files selected. Exiting.")
        input("Press Enter to close this window...") 
        sys.exit(1)
        
    print(f"Found {len(files_to_process)} file(s) to process. Starting batch update...")
    
    for filepath in files_to_process:
        process_tremer_file(filepath)
        
    print("Batch processing complete!")
