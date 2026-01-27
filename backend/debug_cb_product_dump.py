import urllib.request
import ssl

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

url = "https://www.creative-biolabs.com/bsab/bispecific-ang2-vegfa-tandem-diabody-12729.htm"
headers = {'User-Agent': 'Mozilla/5.0'}

req = urllib.request.Request(url, headers=headers)
try:
    with urllib.request.urlopen(req, context=ctx) as response:
        content = response.read().decode('utf-8')
        with open("backend/cb_product_dump_python.html", "w", encoding="utf-8") as f:
            f.write(content)
    print("Success")
except Exception as e:
    print(f"Error: {e}")
