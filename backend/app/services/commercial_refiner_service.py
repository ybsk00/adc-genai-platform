import asyncio
import logging
from typing import List, Dict, Any

from app.core.supabase import supabase
from app.services.ambeed_crawler import ambeed_crawler
from app.services.rag_service import rag_service

logger = logging.getLogger(__name__)

async def batch_refine_commercial_reagents(limit: int = 50, mode: str = "partial"):
    """
    Batch refine 'commercial_reagents' table.
    Target: Records where ai_refined is False (or null).
    Action:
      1. Fetch data.
      2. Call AmbeedCrawler._enrich_with_gemini (or logic).
      3. Generate Embedding (768d).
      4. Update DB.
    """
    logger.info(f"🚀 [Commercial Batch] Starting refinement (Limit: {limit}, Mode: {mode})...")
    
    try:
        # 1. Fetch unrefined records
        # Note: 'ai_refined' column might not exist or be false.
        query = supabase.table("commercial_reagents")\
            .select("*")\
            .eq("ai_refined", False)
            
        if mode == "partial":
            query = query.limit(limit)
        
        response = query.execute()
        items = response.data
        
        if not items:
            logger.info("✨ No pending commercial reagents found.")
            return {"count": 0, "message": "No pending items"}
            
        logger.info(f"🔍 Found {len(items)} items to refine.")
        
        processed_count = 0
        error_count = 0
        
        # 2. Loop & Process
        for item in items:
            try:
                # Prepare text for AI
                description = item.get("body_text") or item.get("summary") or item.get("product_name")
                
                # Ensure properties is a dict
                current_props = item.get("properties") or {}
                item["properties"] = current_props

                # Call AI Enrichment
                ai_data = await ambeed_crawler._enrich_with_gemini(description, item)
                
                if not ai_data:
                    ai_data = {}

                # Merge AI data
                target = ai_data.get("target")
                summary = ai_data.get("summary")
                ai_props = ai_data.get("properties") or {}
                
                # Update local item dict for embedding generation
                item["target"] = target
                item["summary"] = summary
                
                # Generate Embedding
                embedding_text = f"{item.get('product_name', '')} {target or ''} {summary or ''} {item.get('smiles_code', '')}"
                embedding = await rag_service.generate_embedding(embedding_text)
                
                # Prepare DB Update
                update_data = {
                    "target": target,
                    "summary": summary,
                    "ai_refined": True,
                    "embedding": embedding,
                    "properties": {**current_props, **ai_props, "ai_analysis": ai_data}
                }
                
                # Update DB
                supabase.table("commercial_reagents").update(update_data).eq("id", item["id"]).execute()
                processed_count += 1
                logger.info(f"   ✅ Refined: {item.get('product_name', 'Unknown')[:40]}...")
                
                # Rate limit safety
                await asyncio.sleep(0.5)
                
            except Exception as e:
                logger.error(f"   ❌ Error processing {item.get('id')}: {e}")
                error_count += 1
        
        logger.info(f"🎉 [Commercial Batch] Completed. Processed: {processed_count}, Errors: {error_count}")
        return {"count": processed_count, "errors": error_count}

    except Exception as e:
        logger.error(f"❌ Batch Refine Failed: {e}")
        return {"error": str(e)}