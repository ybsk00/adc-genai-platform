import os
import google.generativeai as genai
from dotenv import load_dotenv

load_dotenv("backend/.env")
api_key = os.getenv("GOOGLE_API_KEY")
genai.configure(api_key=api_key)

print("Checking available models...")
try:
    for m in genai.list_models():
        if 'generateContent' in m.supported_generation_methods:
            print(m.name)
except Exception as e:
    print(f"Error listing models: {e}")

print("\nTesting gemini-2.0-flash...")
try:
    model = genai.GenerativeModel('gemini-2.0-flash')
    response = model.generate_content("test")
    print("SUCCESS: gemini-2.0-flash exists and works.")
except Exception as e:
    print(f"FAILED: gemini-2.0-flash failed: {e}")
