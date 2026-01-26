import asyncio
import logging
import sys
import os

# 현재 경로 설정
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 로그가 화면에 잘 보이도록 설정
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)

from app.services.pubmed_knowledge_service import pubmed_knowledge_service

async def main():
    print("\n" + "="*50)
    print("🚀 배치 프로세스 시작 (목표: 약물 10개)")
    print("="*50)
    
    # 배치 실행
    result = await pubmed_knowledge_service.run_batch(batch_size=10, mode='incremental')
    
    print("\n" + "="*50)
    print("✅ 배치 작업 완료!")
    print(f"   - 상태: {result.get('status')}")
    print(f"   - 처리된 약물 수: {result.get('total_drugs', 0)}")
    print(f"   - 저장된 논문 수: {result.get('papers_saved', 0)}")
    print(f"   - 에러 발생: {result.get('errors', 0)}")
    print("="*50)

if __name__ == "__main__":
    # 윈도우 환경에서 asyncio 루프 정책 설정 (필요시)
    if sys.platform == 'win32':
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    asyncio.run(main())
