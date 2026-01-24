import re
import os

file_path = 'ambeed_debug_search.html'
if not os.path.exists(file_path):
    print(f"File not found: {file_path}")
    exit(1)

with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()
    print(f"File length: {len(content)}")
    print(f"Count of '/products/': {content.count('/products/')}")
    print(f"Count of '/record/': {content.count('/record/')}")
    
    # Find all links containing /products/
    links = re.findall(r'href="([^"]*/products/[^"]*)"', content)
    print(f"Found {len(links)} product links (regex):")
    for link in links[:5]:
        print(f" - {link}")

    # Find all links containing /record/
    links_record = re.findall(r'href="([^"]*/record/[^"]*)"', content)
    print(f"Found {len(links_record)} record links (regex):")
    for link in links_record[:5]:
        print(f" - {link}")
