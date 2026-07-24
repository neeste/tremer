import csv
import glob
import os

all_markers = []
seen = set()

# Process all CSVs in ../YSTR/ (all subdirectories)
csv_files = glob.glob('../YSTR/*/*.csv')
if not csv_files:
    csv_files = glob.glob('../YSTR/*.csv')

multi_copy_counts = {
    'DYS385': 2,
    'DYS459': 2,
    'YCAII': 2,
    'CDY': 2,
    'DYF395S1': 2,
    'DYS413': 2,
    'DYS464': 4,
    'DYF406S1': 1, # Wait, is this multi-copy?
    'DYS389I': 1,
    'DYS389II': 1
}

for f in csv_files:
    with open(f, 'r') as file:
        reader = csv.reader(file)
        try:
            headers = next(reader)
        except StopIteration:
            continue
        for h in headers:
            h = h.strip()
            if not h: continue
            h_lower = h.lower()
            if 'kit' in h_lower or 'name' in h_lower or 'haplogroup' in h_lower or 'snp' in h_lower:
                continue
            
            if h not in seen:
                seen.add(h)
                if h in multi_copy_counts:
                    count = multi_copy_counts[h]
                    for idx in range(count):
                        suffix = chr(ord('a') + idx)
                        all_markers.append(f"{h}{suffix}")
                else:
                    all_markers.append(h)

with open('fallback_labels.h', 'w') as out:
    out.write('#ifndef FALLBACK_LABELS_H\n')
    out.write('#define FALLBACK_LABELS_H\n\n')
    out.write(f'int n_lab = {len(all_markers)};\n')
    out.write(f'const char *str_lab[{len(all_markers)}] = {{\n')
    
    for i in range(0, len(all_markers), 10):
        chunk = all_markers[i:i+10]
        out.write('    ' + ', '.join([f'"{m}"' for m in chunk]) + ',\n')
        
    out.write('};\n\n')
    out.write('#endif\n')

print(f"Generated fallback_labels.h with {len(all_markers)} markers.")
