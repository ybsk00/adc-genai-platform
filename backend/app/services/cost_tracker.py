"""
LLM Cost Tracker Service
LLM 호출 비용 추적 및 일일 한도 관리
"""
import logging
from datetime import datetime, date
from typing import Optional

from app.core.supabase import supabase

logger = logging.getLogger(__name__)


class CostTracker:
    """LLM 비용 추적 및 일일 한도 관리"""
    
    # 가격표 (USD per 1K tokens)
    PRICE_PER_1K_TOKENS = {
        # OpenAI
        "gpt-4o": {"input": 0.005, "output": 0.015},
        "gpt-4o-mini": {"input": 0.00015, "output": 0.0006},
        "gpt-4-turbo": {"input": 0.01, "output": 0.03},
        # Google (Gemini)
        "gemini-2.0-flash": {"input": 0.00001, "output": 0.00004},
        "gemini-2.5-flash": {"input": 0.000075, "output": 0.0003},
        "gemini-1.5-pro": {"input": 0.00125, "output": 0.005},
    }
    
    # 일일 한도 (USD)
    DAILY_LIMIT_USD = 10.0
    
    async def track_usage(
        self, 
        model: str, 
        input_tokens: int, 
        output_tokens: int
    ) -> float:
        """
        사용량 기록 및 비용 계산
        반환: 이번 호출 비용 (USD)
        """
        try:
            # 모델명 정규화 (접두사 제거)
            model_key = model.lower()
            for key in self.PRICE_PER_1K_TOKENS:
                if key in model_key:
                    model_key = key
                    break
            
            prices = self.PRICE_PER_1K_TOKENS.get(model_key, {"input": 0.001, "output": 0.002})
            
            # 비용 계산
            input_cost = (input_tokens / 1000) * prices["input"]
            output_cost = (output_tokens / 1000) * prices["output"]
            total_cost = input_cost + output_cost
            
            # DB에 기록
            supabase.table("llm_usage_logs").insert({
                "model": model,
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "cost_usd": total_cost
            }).execute()
            
            logger.info(f"💰 LLM Cost: ${total_cost:.6f} ({model}, {input_tokens}+{output_tokens} tokens)")
            
            return total_cost
            
        except Exception as e:
            logger.error(f"Cost tracking error: {e}")
            return 0.0
    
    async def get_daily_usage(self) -> float:
        """오늘 사용한 총 비용 조회 (USD)"""
        try:
            today = date.today().isoformat()
            
            result = supabase.table("llm_usage_logs")\
                .select("cost_usd")\
                .gte("created_at", f"{today}T00:00:00")\
                .execute()
            
            total = sum(row["cost_usd"] for row in result.data) if result.data else 0.0
            return total
            
        except Exception as e:
            logger.error(f"Get daily usage error: {e}")
            return 0.0
    
    async def is_over_limit(self) -> bool:
        """일일 한도 초과 여부"""
        daily_usage = await self.get_daily_usage()
        is_over = daily_usage >= self.DAILY_LIMIT_USD
        
        if is_over:
            logger.warning(f"⚠️ Daily LLM cost limit reached: ${daily_usage:.2f} >= ${self.DAILY_LIMIT_USD}")
        
        return is_over
    
    async def get_usage_summary(self) -> dict:
        """사용량 요약"""
        daily_usage = await self.get_daily_usage()
        return {
            "daily_usage_usd": daily_usage,
            "daily_limit_usd": self.DAILY_LIMIT_USD,
            "remaining_usd": max(0, self.DAILY_LIMIT_USD - daily_usage),
            "is_over_limit": daily_usage >= self.DAILY_LIMIT_USD
        }


# 싱글톤 인스턴스
cost_tracker = CostTracker()
