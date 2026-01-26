from typing import List, Dict, Any
from pydantic import BaseModel
import google.generativeai as genai
import asyncio
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

    async def search_knowledge_base(self, query: str, top_k: int = 3) -> List[Dict[str, Any]]:
        """Search Knowledge Base (Gemini 768-dim)"""
        try:
            embedding = await self.generate_embedding(query)
            response = supabase.rpc(
                "match_knowledge_base",
                {
                    "query_embedding": embedding,
                    "match_threshold": 0.5,
                    "match_count": top_k
                }
            ).execute()
            # Add source type metadata
            results = []
            for item in response.data:
                item['source'] = 'knowledge_base'
                results.append(item)
            return results
        except Exception as e:
            print(f"Error searching Knowledge Base: {str(e)}")
            return []

    async def search_antibody_library(self, query: str, top_k: int = 3) -> List[Dict[str, Any]]:
        """Search Antibody Library (Gemini 768-dim)"""
        try:
            embedding = await self.generate_embedding(query)
            response = supabase.rpc(
                "match_antibody_library",
                {
                    "query_embedding": embedding,
                    "match_threshold": 0.5,
                    "match_count": top_k
                }
            ).execute()
            for item in response.data:
                item['source'] = 'antibody_library'
            return response.data
        except Exception as e:
            print(f"Error searching Antibody Library: {str(e)}")
            return []

    async def search_golden_set(self, query: str, top_k: int = 3) -> List[Dict[str, Any]]:
        """Search Golden Set (Gemini 768-dim assumed)"""
        try:
            embedding = await self.generate_embedding(query)
            # Note: Ensure match_golden_set accepts 768-dim vector
            response = supabase.rpc(
                "match_golden_set",
                {
                    "query_embedding": embedding,
                    "match_threshold": 0.5,
                    "match_count": top_k
                }
            ).execute()
            for item in response.data:
                item['source'] = 'golden_set'
            return response.data
        except Exception as e:
            print(f"Error searching Golden Set: {str(e)}")
            return []

    async def search_all(self, query: str, top_k_per_source: int = 3) -> List[Dict[str, Any]]:
        """
        Multi-Source Search: Queries all 3 DBs and aggregates results.
        Returns a unified list of relevant chunks.
        """
        # Run searches in parallel
        results = await asyncio.gather(
            self.search_golden_set(query, top_k_per_source),
            self.search_knowledge_base(query, top_k_per_source),
            self.search_antibody_library(query, top_k_per_source)
        )
        
        # Flatten and sort by similarity
        all_results = []
        for res in results:
            all_results.extend(res)
            
        # Sort by similarity desc
        all_results.sort(key=lambda x: x.get('similarity', 0), reverse=True)
        
        return all_results

rag_service = RAGService()
