import asyncio
import os
import sys
import logging

# Add backend directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), 'app'))

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from app.core.supabase import supabase
from app.services.ai_refiner import ai_refiner

async def debug_refiner():
    with open("debug_refiner_output.txt", "w", encoding="utf-8") as f:
        f.write("--- Debugging AI Refiner ---\n")
        
        # 1. Fetch a sample record
        try:
            # Check if column exists first by selecting it
            try:
                f.write("Checking if ai_refined column exists...\n")
                supabase.table("commercial_reagents").select("ai_refined").limit(1).execute()
                f.write("Column ai_refined exists.\n")
            except Exception as e:
                f.write(f"Column check failed: {e}\n")

            f.write("Fetching pending record...\n")
            # Try to find records where ai_refined is False (default) or null
            res = supabase.table("commercial_reagents").select("*").eq("ai_refined", False).limit(1).execute()
            if not res.data:
                f.write("No pending records (ai_refined=False) found. Trying null...\n")
                res = supabase.table("commercial_reagents").select("*").is_("ai_refined", "null").limit(1).execute()
            
            if not res.data:
                f.write("No pending records found in commercial_reagents.\n")
                return
                
            record = res.data[0]
            f.write(f"Testing with Record ID: {record['id']}\n")
            f.write(f"Product Name: {record['product_name']}\n")
            
            # 2. Run Refiner
            f.write("Running refine_single_record...\n")
            result = await ai_refiner.refine_single_record(record)
            
            f.write("\n[Result]\n")
            f.write(str(result) + "\n")
            
        except Exception as e:
            f.write(f"\n[Error] {e}\n")
            import traceback
            traceback.print_exc(file=f)

if __name__ == "__main__":
    asyncio.run(debug_refiner())
