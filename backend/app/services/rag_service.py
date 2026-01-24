from typing import List, Dict, Any
from pydantic import BaseModel
import google.generativeai as genai
from app.core.config import settings
from app.core.supabase import supabase

class RAGService:
    def __init__(self):
        self.api_key = settings.GOOGLE_API_KEY
        self.model_id = "models/text-embedding-004"
        if self.api_key:
            genai.configure(api_key=self.api_key)
        else:
            print("Warning: GOOGLE_API_KEY is missing. RAG Service will not function.")

    async def generate_embedding(self, text: str) -> List[float]:
        """Generate embedding for text using Gemini (768 dimensions)"""
        if not self.api_key:
            print("Error: Gemini API key not available")
            return []
            
        try:
            # text-embedding-004 supports 'output_dimensionality'
            result = genai.embed_content(
                model=self.model_id,
                content=text,
                task_type="retrieval_document",
                output_dimensionality=768
            )
            return result['embedding']
        except Exception as e:
            print(f"Error generating Gemini embedding: {str(e)}")
            return []


    async def index_golden_set_item(self, item_id: str, description: str, properties: Dict[str, Any]):
        """
        Index a Golden Set item into the vector database.
        """
        try:
            # 1. Generate Embedding
            # Combine description and key properties for better retrieval
            content = f"{description}\n"
            for k, v in properties.items():
                if v and isinstance(v, (str, int, float)):
                    content += f"{k}: {v}\n"
            
            embedding = await self.generate_embedding(content)
            
            # 2. Store in Supabase
            data = {
                "source_id": item_id,
                "chunk_content": content,
                "embedding": embedding
            }
            
            supabase.table("golden_set_embeddings").insert(data).execute()
            print(f"Successfully indexed golden set item: {item_id}")
            
        except Exception as e:
            print(f"Error indexing golden set item {item_id}: {str(e)}")
            # TODO: Log error to DB

    async def search(self, query: str, top_k: int = 5) -> List[Dict[str, Any]]:
        """
        Search the vector database for relevant chunks.
        """
        try:
            embedding = await self.generate_embedding(query)
            
            # Call the RPC function defined in init.sql
            response = supabase.rpc(
                "match_golden_set",
                {
                    "query_embedding": embedding,
                    "match_threshold": 0.5,
                    "match_count": top_k
                }
            ).execute()
            
            return response.data
            
        except Exception as e:
            print(f"Error searching vector DB: {str(e)}")
            return []

rag_service = RAGService()
