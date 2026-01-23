"""Quick verification test for PubMed Knowledge Service"""
import asyncio
import sys
sys.path.insert(0, '.')

async def main():
    from app.services.pubmed_knowledge_service import PubMedKnowledgeService
    from app.core.supabase import supabase
    
    service = PubMedKnowledgeService()
    
    # Test 1: Get drugs
    print("=" * 50)
    print("Test 1: Get target drugs")
    drugs = await service.get_target_drugs(limit=5)
    print(f"Found {len(drugs)} drugs")
    if drugs:
        for d in drugs[:3]:
            print(f"  - {d['name']}")
    
    # Test 2: Search PubMed  
    print("\n" + "=" * 50)
    print("Test 2: Search PubMed for 'Enhertu'")
    articles = await service.search_pubmed_for_drug("Enhertu", max_results=2)
    print(f"Found {len(articles)} articles")
    if articles:
        print(f"  First article: {articles[0]['title'][:60]}...")
    
    # Test 3: Check KB
    print("\n" + "=" * 50)
    print("Test 3: Knowledge Base stats")
    res = supabase.table("knowledge_base").select("count", count="exact").eq("source_type", "PubMed").execute()
    print(f"Current PubMed records: {res.count}")
    
    # Test 4: Mini batch (1 drug only)
    print("\n" + "=" * 50)
    print("Test 4: Process single drug")
    if drugs:
        saved = await service.process_single_drug(drugs[0])
        print(f"Saved {saved} articles for {drugs[0]['name']}")
    
    print("\n" + "=" * 50)
    print("DONE!")

if __name__ == "__main__":
    asyncio.run(main())
