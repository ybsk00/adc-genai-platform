import os
from typing import List, Dict, Any
from pydantic import BaseModel
from openai import AsyncOpenAI
from app.core.supabase import supabase

class RAGService:
    def __init__(self):
        self.api_key = os.getenv("OPENAI_API_KEY")
        self.client = None

    def _get_client(self):
        if not self.client:
            if not self.api_key:
                print("Warning: OPENAI_API_KEY is missing. RAG Service will not function.")
                return None
            self.client = AsyncOpenAI(api_key=self.api_key)
        return self.client

    async def generate_embedding(self, text: str) -> List[float]:
        """Generate embedding for text using OpenAI"""
        client = self._get_client()
        if not client:
            print("Error: OpenAI client not available")
            return []
            
        response = await client.embeddings.create(
            input=text,
            model="text-embedding-3-small"
        )
        return response.data[0].embedding


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
