import os
import google.generativeai as genai
from dotenv import load_dotenv

load_dotenv("backend/.env")

api_key = os.getenv("GOOGLE_API_KEY")
if not api_key:
    print("GOOGLE_API_KEY not found")
    exit(1)

genai.configure(api_key=api_key)

# Test with gemini-2.0-flash
model = genai.GenerativeModel('gemini-2.0-flash')

prompt = "Hello, tell me about ADC drugs in one sentence."

try:
    print(f"Testing Gemini 2.0 Flash with prompt: {prompt}")
    response = model.generate_content(prompt)
    print(f"RESPONSE: {response.text}")
except Exception as e:
    print(f"ERROR: {e}")

# Test with structured prompt
structured_prompt = """You are a Clinical Trial Analyst.
Analyze this: Title: A Phase 1 Study of ADC.
Output JSON: {"drug_name": "ADC"}"""

try:
    print(f"\nTesting with structured prompt: {structured_prompt}")
    response = model.generate_content(structured_prompt)
    print(f"RESPONSE: {response.text}")
except Exception as e:
    print(f"ERROR: {e}")
