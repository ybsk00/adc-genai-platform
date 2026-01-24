import os
import asyncio
import logging
import json
from typing import List, Dict, Any
from dotenv import load_dotenv

# 환경 변수 로드
env_path = os.path.join(os.path.dirname(__file__), '.env')
load_dotenv(env_path)

from app.core.config import settings
from app.core.supabase import supabase
from app.services.rag_service import rag_service

# 로깅 설정
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class Reindexer:
    def __init__(self):
        self.batch_size = 20
        self.delay_between_batches = 2.0  # API 할당량 고려 (TPM/RPM)

    async def reindex_golden_set(self):
        logger.info("🚀 [Golden Set] 재인덱싱 시작...")
        try:
            # 1. 기존 임베딩 삭제 (768차원으로 새로 채우기 위해)
            supabase.table("golden_set_embeddings").delete().neq("id", "00000000-0000-0000-0000-000000000000").execute()
            
            # 2. 원본 데이터 가져오기
            response = supabase.table("golden_set_library").select("*").eq("status", "approved").execute()
            items = response.data
            logger.info(f"Found {len(items)} approved golden set items.")

            for i in range(0, len(items), self.batch_size):
                batch = items[i:i + self.batch_size]
                for item in batch:
                    props = item.get("properties", {}) or {}
                    await rag_service.index_golden_set_item(
                        item_id=item["id"],
                        description=item.get("description", ""),
                        properties=props
                    )
                logger.info(f"Indexed {min(i + self.batch_size, len(items))}/{len(items)} golden set items.")
                await asyncio.sleep(self.delay_between_batches)
            
            logger.info("✅ [Golden Set] 재인덱싱 완료.")
        except Exception as e:
            logger.error(f"❌ Golden Set 재인덱싱 실패: {e}")

    async def reindex_commercial_reagents(self):
        logger.info("🚀 [Commercial Reagents] 재인덱싱 시작...")
        try:
            response = supabase.table("commercial_reagents").select("id, product_name, target, properties, summary").execute()
            items = response.data
            logger.info(f"Found {len(items)} commercial reagents.")

            for i in range(0, len(items), self.batch_size):
                batch = items[i:i + self.batch_size]
                for item in batch:
                    # 임베딩용 텍스트 결합
                    content = f"Product: {item['product_name']}\nTarget: {item.get('target', 'N/A')}\n"
                    if item.get("summary"):
                        content += f"Summary: {item['summary']}\n"
                    props = item.get("properties", {}) or {}
                    for k, v in props.items():
                        if v: content += f"{k}: {v}\n"
                    
                    embedding = await rag_service.generate_embedding(content)
                    if embedding:
                        supabase.table("commercial_reagents").update({"embedding": embedding}).eq("id", item["id"]).execute()
                
                logger.info(f"Indexed {min(i + self.batch_size, len(items))}/{len(items)} commercial reagents.")
                await asyncio.sleep(self.delay_between_batches)
            
            logger.info("✅ [Commercial Reagents] 재인덱싱 완료.")
        except Exception as e:
            logger.error(f"❌ Commercial Reagents 재인덱싱 실패: {e}")

    async def reindex_knowledge_base(self):
        logger.info("🚀 [Knowledge Base] 재인덱싱 시작...")
        try:
            # knowledge_base 테이블은 content(초록)를 기반으로 임베딩
            response = supabase.table("knowledge_base").select("id, title, content, abstract").execute()
            items = response.data
            logger.info(f"Found {len(items)} knowledge base items.")

            for i in range(0, len(items), self.batch_size):
                batch = items[i:i + self.batch_size]
                for item in batch:
                    # Use abstract if available, otherwise content
                    text_content = item.get('abstract') or item.get('content') or ""
                    content = f"Title: {item['title']}\nContent: {text_content}"
                    embedding = await rag_service.generate_embedding(content)
                    if embedding:
                        # knowledge_base 테이블에 embedding 컬럼이 있는지 확인 필요 (없으면 추가 로직 필요할 수 있음)
                        try:
                            supabase.table("knowledge_base").update({"embedding": embedding}).eq("id", item["id"]).execute()
                        except Exception as inner_e:
                            logger.warning(f"Skipping KB item {item['id']} (Check if embedding column exists): {inner_e}")
                
                logger.info(f"Indexed {min(i + self.batch_size, len(items))}/{len(items)} knowledge base items.")
                await asyncio.sleep(self.delay_between_batches)
            
            logger.info("✅ [Knowledge Base] 재인덱싱 완료.")
        except Exception as e:
            logger.error(f"❌ Knowledge Base 재인덱싱 실패: {e}")

    async def run_all(self):
        await self.reindex_golden_set()
        await self.reindex_commercial_reagents()
        await self.reindex_knowledge_base()

if __name__ == "__main__":
    reindexer = Reindexer()
    asyncio.run(reindexer.run_all())
