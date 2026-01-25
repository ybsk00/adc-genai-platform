import asyncio
import json
import os
import google.generativeai as genai
from dotenv import load_dotenv

# Load env from backend/.env
load_dotenv("backend/.env")

class Settings:
    GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")

settings = Settings()

async def test_chat():
    try:
        genai.configure(api_key=settings.GOOGLE_API_KEY)
        
        # Try the same logic as in admin.py
        try:
            print("Trying gemini-2.0-flash...")
            model = genai.GenerativeModel('gemini-2.0-flash')
            # Test if it actually works
            response = model.generate_content("test")
            print("gemini-2.0-flash worked")
        except Exception as e:
            print(f"⚠️ gemini-2.0-flash failed: {e}")
            print("Falling back to gemini-1.5-flash...")
            model = genai.GenerativeModel('gemini-1.5-flash')
        
        context_str = json.dumps({"test": "data"}, indent=2, ensure_ascii=False)
        
        system_prompt = f"You are an assistant. RAW DATA: {context_str}"
        full_prompt = f"{system_prompt}\n\nUser Question: hello"
        
        print(f"Calling generate_content with model: {model.model_name}")
        response = await asyncio.get_event_loop().run_in_executor(
            None, lambda: model.generate_content(full_prompt)
        )
        
        print(f"ANSWER: {response.text.strip()}")
        
    except Exception as e:
        print(f"AI Chat Error: {e}")

if __name__ == "__main__":
    asyncio.run(test_chat())
