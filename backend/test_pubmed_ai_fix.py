"""
Verification Script for PubMed AI Fix
Tests:
1. Connectivity (Single API call)
2. Mylotarg Sample (Target extraction verification)
3. Batch Log Watch (10 items)
"""
import os
import sys
import asyncio
import logging
from dotenv import load_dotenv

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

load_dotenv()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("verification_logs.txt", mode='w', encoding='utf-8')
    ]
)
logger = logging.getLogger(__name__)

# Redirect print to logger
def print(*args, **kwargs):
    msg = " ".join(map(str, args))
    logger.info(msg)

from app.services.pubmed_knowledge_service import pubmed_knowledge_service

async def test_step_1_connectivity():
    print("\n🧪 [Step 1] Testing Connectivity (Single API Call)...")
    dummy_abstract = """
    Trastuzumab deruxtecan (T-DXd) is an antibody-drug conjugate targeting HER2. 
    In this Phase 3 trial for breast cancer, T-DXd showed significant improvement in progression-free survival compared to T-DM1.
    The safety profile was consistent with previous studies.
    """
    dummy_title = "Trastuzumab deruxtecan in HER2-positive breast cancer"
    
    result = await pubmed_knowledge_service.analyze_with_gemini(dummy_abstract, dummy_title, "Trastuzumab deruxtecan")
    
    print(f"   Result: {result}")
    
    if result.get("target") == "HER2" and result.get("relevance_score", 0) > 0.8:
        print("✅ Step 1 Passed: 404 Error Gone, Analysis Successful.")
        return True
    else:
        print("❌ Step 1 Failed: Invalid analysis result.")
        return False

async def test_step_2_mylotarg():
    print("\n🧪 [Step 2] Testing Mylotarg (Gemtuzumab ozogamicin)...")
    
    # Clean up existing Mylotarg data to avoid duplicate skipping
    from app.core.supabase import supabase
    print("   🧹 Cleaning up existing Mylotarg papers...")
    supabase.table("knowledge_base").delete().ilike("title", "%Gemtuzumab%").execute()
    
    # Mock drug data
    mylotarg = {
        "id": "test_mylotarg",
        "name": "Gemtuzumab ozogamicin",
        "generic_name": "Gemtuzumab ozogamicin"
    }
    
    # Process single drug with debug logging
    print("   🔍 Searching PubMed for Gemtuzumab ozogamicin...")
    articles = await pubmed_knowledge_service.search_pubmed_for_drug("Gemtuzumab ozogamicin", max_results=5)
    print(f"   📚 Found {len(articles)} articles.")
    
    if not articles:
        print("   ❌ No articles found. Check search query or date filter.")
        return False
    
    # Clean up specific titles to ensure they are re-processed
    from app.core.supabase import supabase
    print("   🧹 Cleaning up specific titles found...")
    for article in articles:
        supabase.table("knowledge_base").delete().eq("title", article["title"]).execute()
        
    saved_count = 0
    for article in articles:
        print(f"   📄 Processing: {article['title'][:50]}...")
        # Check duplicate
        if await pubmed_knowledge_service.is_duplicate(article["title"]):
            print("      ⚠️ Duplicate (Skipped)")
            continue
            
        analysis = await pubmed_knowledge_service.analyze_with_gemini(article["abstract"], article["title"], "Gemtuzumab ozogamicin")
        if await pubmed_knowledge_service.save_to_knowledge_base(article, analysis, mylotarg["id"]):
            saved_count += 1
            print(f"      ✅ Saved (Score: {analysis.get('relevance_score')})")
    
    print(f"   Saved {saved_count} papers for Mylotarg.")
    
    if saved_count > 0:
        print("✅ Step 2 Passed: Mylotarg processed successfully.")
        return True
    else:
        print("❌ Step 2 Failed: No papers saved for Mylotarg.")
        return False

async def test_step_3_batch_log_watch():
    print("\n🧪 [Step 3] Log Watch (Batch of 10)...")
    
    # Run a small batch
    result = await pubmed_knowledge_service.run_batch(batch_size=10, mode="incremental")
    
    print(f"   Batch Result: {result}")
    
    if result["status"] == "completed" and result["errors"] == 0:
        print("✅ Step 3 Passed: Batch completed with 0 errors.")
        return True
    else:
        print(f"❌ Step 3 Failed: Errors detected ({result.get('errors')}).")
        return False

async def main():
    print("🚀 Starting 3-Step Verification...")
    
    if not await test_step_1_connectivity():
        print("⛔ Aborting after Step 1 failure.")
        return
        
    if not await test_step_2_mylotarg():
        print("⛔ Aborting after Step 2 failure.")
        return
        
    await test_step_3_batch_log_watch()
    
    print("\n🎉 All Tests Completed.")

if __name__ == "__main__":
    asyncio.run(main())
