import re
import os

file_path = 'ambeed_debug_search.html'
if not os.path.exists(file_path):
    print(f"File not found: {file_path}")
    exit(1)

with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()
    
    matches = [m.start() for m in re.finditer('ADC', content)]
    with open('adc_matches.txt', 'w', encoding='utf-8') as out_f:
        out_f.write(f"Found {len(matches)} matches for 'ADC'\n")
        
        for idx in matches:
            start = max(0, idx - 100)
            end = min(len(content), idx + 100)
            context = content[start:end]
            out_f.write(f"Context at {idx}: {context}\n")
            out_f.write("-" * 20 + "\n")
    print("Matches written to adc_matches.txt")
