"""
OpenFDA Targeted Search 테스트 스크립트
- ClinicalTrials 약물명을 사용한 정밀 검색 테스트
"""
import asyncio
import sys
sys.path.insert(0, '.')

from app.services.openfda_service import openfda_service

async def main():
    print("=" * 60)
    print("🧪 OpenFDA Targeted Search Test")
    print("=" * 60)
    
    # Test 1: 기본 쿼리만 (Targeted Search 없이)
    print("\n📊 Test 1: 기본 쿼리만 (Targeted Search OFF)")
    print("-" * 40)
    results_basic = await openfda_service.fetch_all_approved_adcs(
        limit=100,
        mode="full",
        use_broad_search=True,
        use_targeted_search=False
    )
    print(f"✅ Basic Search Results: {len(results_basic)} records")
    
    # Test 2: Targeted Search 포함 (ClinicalTrials 약물명 사용)
    print("\n📊 Test 2: Targeted Search ON (ClinicalTrials 약물명)")
    print("-" * 40)
    results_targeted = await openfda_service.fetch_all_approved_adcs(
        limit=100,
        mode="full",
        use_broad_search=True,
        use_targeted_search=True
    )
    print(f"✅ Targeted Search Results: {len(results_targeted)} records")
    
    # Summary
    print("\n" + "=" * 60)
    print("📈 SUMMARY")
    print("=" * 60)
    print(f"  Basic Search:    {len(results_basic)} records")
    print(f"  Targeted Search: {len(results_targeted)} records")
    improvement = len(results_targeted) - len(results_basic)
    pct = ((len(results_targeted) / max(len(results_basic), 1)) - 1) * 100
    print(f"  Improvement:     +{improvement} records ({pct:.1f}% increase)")
    
    # Sample data preview
    if results_targeted:
        print("\n📋 Sample Records (first 10):")
        for i, r in enumerate(results_targeted[:10]):
            openfda = r.get("openfda", {})
            brand = openfda.get("brand_name", ["N/A"])[0] if openfda.get("brand_name") else "N/A"
            generic = openfda.get("generic_name", ["N/A"])[0] if openfda.get("generic_name") else "N/A"
            print(f"  {i+1}. {brand} ({generic})")

if __name__ == "__main__":
    asyncio.run(main())
