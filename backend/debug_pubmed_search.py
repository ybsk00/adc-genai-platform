import asyncio
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from app.services.pubmed_knowledge_service import pubmed_knowledge_service

async def debug_drugs():
    print("🔍 검사 대상 약물 확인...")
    drugs = await pubmed_knowledge_service.get_target_drugs(limit=10)
    
    for i, drug in enumerate(drugs):
        name = drug['name']
        print(f"[{i+1}] {name}")
        
        # 실제 검색 쿼리 테스트
        query = pubmed_knowledge_service.build_search_query(name)
        print(f"   - Query: {query[:100]}...")
        
        # 검색 결과 개수 확인
        articles = await pubmed_knowledge_service.search_pubmed_for_drug(name, max_results=1)
        print(f"   - 결과 수: {len(articles)}개")

if __name__ == "__main__":
    if sys.platform == 'win32':
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    asyncio.run(debug_drugs())
