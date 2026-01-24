import os

log_file = 'crawler_debug_v4.log'
if os.path.exists(log_file):
    try:
        with open(log_file, 'r', encoding='utf-16') as f:
            for line in f:
                if "Found: Cat=" in line:
                    print(line.strip())
    except UnicodeError:
        try:
            with open(log_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if "Found: Cat=" in line:
                        print(line.strip())
        except Exception as e:
            print(f"Error reading log: {e}")
else:
    print("Log file not found.")
