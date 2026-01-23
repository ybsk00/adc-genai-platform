"""
PubMed Knowledge Service 테스트 스크립트
단위 테스트 및 샘플 실행
"""
import asyncio
import sys
import os

# 경로 설정
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from app.core.supabase import supabase
from app.services.pubmed_knowledge_service import PubMedKnowledgeService


async def test_get_target_drugs():
    """1. 약물 리스트 추출 테스트"""
    print("\n" + "="*60)
    print("📋 테스트 1: 약물 리스트 추출")
    print("="*60)
    
    service = PubMedKnowledgeService()
    drugs = await service.get_target_drugs(limit=10)
    
    print(f"✅ 추출된 약물 수: {len(drugs)}")
    for drug in drugs[:5]:
        print(f"   - {drug['name']} (Generic: {drug.get('generic_name', 'N/A')})")
    
    return len(drugs) > 0


async def test_build_search_query():
    """2. 검색 쿼리 생성 테스트"""
    print("\n" + "="*60)
    print("🔍 테스트 2: PubMed 검색 쿼리 생성")
    print("="*60)
    
    service = PubMedKnowledgeService()
    
    test_drugs = ["Enhertu", "Trastuzumab deruxtecan", "DS-8201a"]
    for drug in test_drugs:
        query = service.build_search_query(drug)
        print(f"   {drug} → {query[:80]}...")
    
    return True


async def test_pubmed_search():
    """3. PubMed 검색 테스트"""
    print("\n" + "="*60)
    print("🔬 테스트 3: PubMed 검색 (Enhertu)")
    print("="*60)
    
    service = PubMedKnowledgeService()
    articles = await service.search_pubmed_for_drug("Enhertu", max_results=3)
    
    print(f"✅ 검색된 논문 수: {len(articles)}")
    for article in articles:
        print(f"   - PMID: {article.get('pmid')}")
        print(f"     Title: {article.get('title', '')[:60]}...")
        print(f"     Abstract: {len(article.get('abstract', ''))} chars")
    
    return len(articles) > 0


async def test_gemini_analysis():
    """4. Gemini AI 분석 테스트"""
    print("\n" + "="*60)
    print("🤖 테스트 4: Gemini AI 분석")
    print("="*60)
    
    service = PubMedKnowledgeService()
    
    # 샘플 초록
    sample_abstract = """
    Trastuzumab deruxtecan (T-DXd, DS-8201a) is an antibody-drug conjugate (ADC) 
    composed of a humanized anti-HER2 antibody, a cleavable tetrapeptide-based linker, 
    and a potent topoisomerase I inhibitor payload. In the DESTINY-Breast03 trial, 
    T-DXd demonstrated significantly improved progression-free survival (PFS) compared 
    with trastuzumab emtansine (T-DM1) in patients with HER2-positive metastatic breast 
    cancer previously treated with trastuzumab and a taxane. The objective response rate 
    (ORR) was 79.7% for T-DXd versus 34.2% for T-DM1.
    """
    
    analysis = await service.analyze_with_gemini(
        abstract=sample_abstract,
        title="DESTINY-Breast03 Trial Results",
        drug_name="Trastuzumab deruxtecan"
    )
    
    print(f"   Target: {analysis.get('target')}")
    print(f"   Indication: {analysis.get('indication')}")
    print(f"   Relevance Score: {analysis.get('relevance_score')}")
    print(f"   Summary: {analysis.get('summary', '')[:100]}...")
    print(f"   AI Reasoning: {analysis.get('ai_reasoning', '')[:80]}...")
    
    return analysis.get('relevance_score', 0) > 0.5


async def test_save_to_kb():
    """5. Knowledge Base 저장 테스트 (Dry Run)"""
    print("\n" + "="*60)
    print("💾 테스트 5: Knowledge Base 저장 확인")
    print("="*60)
    
    # 기존 PubMed 레코드 수 확인
    result = supabase.table("knowledge_base")\
        .select("count", count="exact")\
        .eq("source_type", "PubMed")\
        .execute()
    
    print(f"   현재 PubMed 레코드 수: {result.count}")
    
    return True


async def test_sample_batch():
    """6. 샘플 배치 실행 (3개 약물만)"""
    print("\n" + "="*60)
    print("🚀 테스트 6: 샘플 배치 실행 (3개 약물)")
    print("="*60)
    
    service = PubMedKnowledgeService()
    result = await service.run_batch(batch_size=3, mode="incremental")
    
    print(f"   상태: {result.get('status')}")
    print(f"   저장된 논문: {result.get('papers_saved', 0)}")
    print(f"   중복 스킵: {result.get('skipped_duplicates', 0)}")
    print(f"   에러: {result.get('errors', 0)}")
    
    return result.get('status') == 'completed'


async def main():
    print("\n" + "="*60)
    print("🧪 PubMed Knowledge Service 테스트 시작")
    print("="*60)
    
    results = []
    
    # 1. 약물 리스트
    results.append(("약물 리스트 추출", await test_get_target_drugs()))
    
    # 2. 쿼리 생성
    results.append(("검색 쿼리 생성", await test_build_search_query()))
    
    # 3. PubMed 검색
    results.append(("PubMed 검색", await test_pubmed_search()))
    
    # 4. Gemini 분석
    results.append(("Gemini AI 분석", await test_gemini_analysis()))
    
    # 5. KB 저장 확인
    results.append(("Knowledge Base 확인", await test_save_to_kb()))
    
    # 선택적: 샘플 배치
    run_sample = input("\n샘플 배치 실행 (3개 약물)? (y/n): ").strip().lower()
    if run_sample == 'y':
        results.append(("샘플 배치 실행", await test_sample_batch()))
    
    # 결과 요약
    print("\n" + "="*60)
    print("📊 테스트 결과 요약")
    print("="*60)
    
    all_passed = True
    for name, passed in results:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"   {status} - {name}")
        if not passed:
            all_passed = False
    
    print("\n" + ("🎉 모든 테스트 통과!" if all_passed else "⚠️ 일부 테스트 실패"))
    
    return all_passed


if __name__ == "__main__":
    asyncio.run(main())
