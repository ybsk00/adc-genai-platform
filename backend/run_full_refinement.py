import asyncio
import logging
import os
import sys
from dotenv import load_dotenv

# Load .env
# Add backend directory to sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

load_dotenv(os.path.join(current_dir, ".env"))

from app.services.ai_refiner import ai_refiner
from batch_refine_commercial import batch_refine_commercial_reagents

# Logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("FULL_REFINER")

async def run_full_batch():
    print("\n" + "="*60)
    print("🚀 STARTING FULL AI REFINEMENT (Golden Set + Commercial)")
    print("="*60 + "\n")
    
    # 1. Golden Set Library (Clinical Trials / FDA / PubMed / etc.)
    print(">>> Phase 1: Refining Golden Set Library (Target: ~2000 items)")
    try:
        # process_pending_records updates DB directly
        await ai_refiner.process_pending_records(max_records=2500, source_filter=None) 
    except Exception as e:
        print(f"❌ Phase 1 Error: {e}")
        
    print("\n" + "-"*60 + "\n")

    # 2. Commercial Reagents (Ambeed / Creative Biolabs)
    print(">>> Phase 2: Refining Commercial Reagents (Target: Any pending)")
    try:
        # custom batch script for commercial table
        res = await batch_refine_commercial_reagents(limit=2500, mode="full")
        print(f"Phase 2 Result: {res}")
    except Exception as e:
        print(f"❌ Phase 2 Error: {e}")

    print("\n" + "="*60)
    print("✅ FULL REFINEMENT JOB COMPLETE")
    print("="*60 + "\n")

if __name__ == "__main__":
    asyncio.run(run_full_batch())
