import argparse
import asyncio
import logging
import sys
import os

# 현재 경로 설정
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 로그 설정
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)

from app.services.pubmed_knowledge_service import pubmed_knowledge_service

async def main():
    parser = argparse.ArgumentParser(description="Golden Set Enricher - PubMed Crawler")
    parser.add_argument("--mode", type=str, choices=["incremental", "discovery"], required=True, help="Mode: incremental (existing drugs) or discovery (new papers)")
    parser.add_argument("--limit", type=int, default=2000, help="Number of papers/drugs to process")
    
    args = parser.parse_args()
    
    print("\n" + "="*60)
    print(f"🚀 Golden Set Enricher 시작")
    print(f"   Mode: {args.mode}")
    print(f"   Limit: {args.limit}")
    print("="*60)
    
    if args.mode == "discovery":
        print("📡 Discovery Mode: 'Mega-Net' Strategy Activated")
        print("   - Yearly Chunking: 2026 ~ 2010 (Reverse Loop)")
        print("   - Expanded Queries: *-mab, *-tin, *-can, NCT extraction")
    else:
        print("💊 Incremental Mode: 기존 약물 데이터 보강...")

    result = await pubmed_knowledge_service.run_batch(batch_size=args.limit, mode=args.mode)
    
    print("\n" + "="*60)
    print("✅ 작업 완료")
    print(f"   Status: {result.get('status')}")
    print(f"   Papers Saved: {result.get('papers_saved', 0)}")
    print(f"   Total Processed: {result.get('processed', 0) if args.mode == 'discovery' else result.get('total_drugs', 0)}")
    if result.get('errors'):
        print(f"   Errors: {result.get('errors')}")
    print("="*60)

if __name__ == "__main__":
    if sys.platform == 'win32':
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    asyncio.run(main())