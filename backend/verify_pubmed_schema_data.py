"""
Verify PubMed Knowledge Base Data Structure
Checks for:
1. Embedding vector presence
2. Structured properties (ORR, Linker, Payload)
3. Gemini Analysis content
"""
import asyncio
import sys
import json
sys.path.insert(0, '.')

from app.core.supabase import supabase

async def verify_latest_pubmed_entry():
    print("🔍 Verifying Knowledge Base Data Structure...")
    
    # Fetch latest PubMed entry
    response = supabase.table("knowledge_base")\
        .select("*")\
        .eq("source_type", "PubMed")\
        .order("created_at", desc=True)\
        .limit(1)\
        .execute()
        
    if not response.data:
        print("❌ No PubMed entries found in knowledge_base table.")
        print("   -> Run the service first to generate data.")
        return
        
    record = response.data[0]
    print(f"✅ Found latest record: {record['title'][:50]}...")
    print(f"   ID: {record['id']}")
    print(f"   Created At: {record['created_at']}")
    
    # 1. Check Embedding
    embedding = record.get("embedding")
    if embedding:
        # Vector might be returned as string or list depending on driver
        if isinstance(embedding, str):
            print(f"✅ Embedding present (String format, length: {len(embedding)})")
        elif isinstance(embedding, list):
            print(f"✅ Embedding present (Vector format, dim: {len(embedding)})")
        else:
            print("✅ Embedding field is present.")
    else:
        print("❌ Embedding is MISSING or NULL.")

    # 2. Check Properties
    props = record.get("properties") or {}
    print("\n📋 Structured Properties Check:")
    keys_to_check = ["orr", "linker", "payload", "outcome_type", "drug_standard_name"]
    
    found_any = False
    for k in keys_to_check:
        val = props.get(k)
        status = f"'{val}'" if val else "Null (Expected if not in text)"
        print(f"   - {k}: {status}")
        if val: found_any = True
        
    if found_any:
        print("✅ Structured extraction is WORKING.")
    else:
        print("⚠️ No specific properties extracted in this sample (might be normal for some abstracts).")
        
    # 3. Check AI Reasoning
    reasoning = record.get("ai_reasoning")
    if reasoning:
        print(f"\n🧠 AI Reasoning: {reasoning[:100]}...")
    else:
        print("\n❌ AI Reasoning is MISSING.")

    print("\n" + "="*50)
    print("Verification Complete.")

if __name__ == "__main__":
    asyncio.run(verify_latest_pubmed_entry())
