"""Debug script to check knowledge_base data and test recovery"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dotenv import load_dotenv
load_dotenv()

from supabase import create_client
import google.generativeai as genai

# Connect to Supabase
client = create_client(
    os.getenv("SUPABASE_URL"),
    os.getenv("SUPABASE_SERVICE_KEY")
)

print("1. Checking knowledge_base for PubMed data with score=0...")
result = client.table("knowledge_base")\
    .select("id, title, content, source_type, relevance_score")\
    .eq("source_type", "PubMed")\
    .eq("relevance_score", 0)\
    .limit(3)\
    .execute()

print(f"   Found {len(result.data)} items")

if result.data:
    for i, item in enumerate(result.data[:3]):
        print(f"\n   Item {i+1}:")
        print(f"   - ID: {item.get('id', 'N/A')}")
        print(f"   - Title: {item.get('title', 'N/A')[:80]}...")
        content = item.get('content', '')
        print(f"   - Content Length: {len(content)} chars")
        if content:
            print(f"   - Content Preview: {content[:200]}...")
else:
    print("   No zero-score PubMed items found!")
    
    # Check if there's any PubMed data at all
    all_pubmed = client.table("knowledge_base")\
        .select("id, relevance_score")\
        .eq("source_type", "PubMed")\
        .limit(10)\
        .execute()
    print(f"\n   Total PubMed items: {len(all_pubmed.data)}")
    if all_pubmed.data:
        scores = [x.get('relevance_score', 0) for x in all_pubmed.data]
        print(f"   Scores: {scores}")

print("\n2. Testing Gemini 2.5 Flash...")
try:
    genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))
    model = genai.GenerativeModel('gemini-2.0-flash-exp')
    
    response = model.generate_content(
        "Return exactly: {\"test\": \"ok\"}",
        safety_settings=[
            {"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}
        ]
    )
    print(f"   Response: {response.text.strip()}")
    print("   ✅ Gemini OK!")
except Exception as e:
    print(f"   ❌ Gemini Error: {e}")

print("\n3. Done!")
